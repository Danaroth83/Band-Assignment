%% MAIN: LOAD DATASETS (REDUCED RESOLUTION - MULTIMODAL)

disp('Reading dataset...');

%% Note: the variable ratio is the desired ratio between the high and low resolution image after processing, not the actual one
if ~exist('ratio','var'), flag_loaddefaultratio=1; else flag_loaddefaultratio=0; end
if ~exist('scale_ref','var'), scale_ref=1; end  % 1=Resizing scale based on MS; 2=based on PAN
if ~exist('PAN_align','var'), PAN_align=2; end  % 1=PAN aligns with MS, 2=PAN aligns with HS
if isequal(PAN_align,'PAN'), PAN_align='HS'; end
if ~exist('flag_imresize_PAN','var'), flag_imresize_PAN=1; end  % 1=Prepare PAN image with imresize, 0=with MTF-filtered resize
if ~exist('flag_imresize_MS','var'), flag_imresize_MS=0; end
if ~exist('flag_imresize_MS_RR','var'), flag_imresize_MS_RR=0; end % 1=Prepare MS image with imresize, 0=with MTF-filtered resize
if ~exist('flag_imresize_PAN_RR','var'), flag_imresize_PAN_RR=1; end
if ~exist('flag_imresize_HS_RR','var'), flag_imresize_HS_RR=0; end
if ~exist('flag_interpolation','var'), flag_interpolation=0; end
if ~exist('flag_avoidEXP','var'), flag_avoidEXP=0; end
if ~exist('flag_avoidPANresize','var'), flag_avoidPANresize=0; end

%% Dataset options
[place,cut_label,Band_label,sensor,sensor_PAN,sensor_MS]=parse_im_tag(im_tag);
ratio_PAN=load_resolution(sensor,im_tag)/load_resolution(sensor_PAN,im_tag,'PAN');
ratio_MS=load_resolution(sensor,im_tag)/load_resolution(sensor_MS,im_tag,'MS');

[Bands_to_sharpen,Bands_to_display]=load_Bands_to_sharpen(sensor,Band_label,place);
[Bands_to_sharpen_MS,Bands_to_display_MS]=load_Bands_to_sharpen(sensor_MS,Band_label,place);
Band_overlap_PAN=load_Band_overlap(sensor,sensor_PAN,'MS','PAN',place);
Band_overlap_MS=load_Band_overlap(sensor,sensor_MS,'MS','MS',place);
[~,L]=load_resolution(sensor,place); % Sensor Radiometric Resolution

if flag_loaddefaultratio==1, ratio=3; end

%%% Cut ROI
[hcut,vcut,edge_cut,~,Qblocks_size]=load_cut(im_tag,cut_label,sensor,sensor_PAN,'PAN');
[~,~,~,edge_cut_MS,~]=load_cut(im_tag,cut_label,sensor,sensor_MS,'MS');
hshift_PAN=0; vshift_PAN=0; hshift_MS=0; vshift_MS=0;

load(['Datasets/',place,'_FR.mat']);
I_PAN=I_MS_LR;
I_MS_LR=I_HS_LR;

[I_PAN,edge_cut_MS]=imcrop_custom(I_PAN,...
    hcut,vcut,ratio_MS,edge_cut_MS+hshift_MS,edge_cut_MS+vshift_MS);

%%% Cut ROI
[I_MS_LR_e,edge_cut]=imcrop_custom(I_MS_LR,hcut,vcut,1,edge_cut,edge_cut);
I_PAN = double(I_PAN);
I_MS_LR_e=double(I_MS_LR_e);

I_MS_loaded=I_MS_LR_e;
I_PAN_loaded=I_PAN;

%% Cutting edges
if ~exist('edge_cut','var'), edge_cut=0; end
if ~exist('ratio_tgt_REF','var'), ratio_tgt_REF=ratio; end

%% Sensor fix
if ~exist('sensor_PAN','var'), sensor_PAN=sensor; end
if ~exist('Bands_to_sharpen','var'), Bands_to_sharpen=1:size(I_MS_loaded,3); end
if ~exist('Bands_to_display','var'), Bands_to_display=1:3; end

%% Manage bands to display
[~,Bands_to_display]=min(abs(repmat(Bands_to_sharpen',[1,3])-repmat(Bands_to_display,[length(Bands_to_sharpen),1])),[],1);
[~,Bands_to_display_MS]=min(abs(repmat(Bands_to_sharpen_MS',[1,3])-repmat(Bands_to_display_MS,[length(Bands_to_sharpen_MS),1])),[],1);

%% Manage bands that overlap
Band_overlap_MS=Band_overlap_MS(Bands_to_sharpen_MS);
% Band_overlap_MS2=unique(find(ismember(Bands_to_sharpen,cell2mat(Band_overlap_MS))));1
Band_overlap_MS=cellfun(@(x) find(ismember(Bands_to_sharpen,x)),Band_overlap_MS,'Un',0);
Band_overlap_PAN=unique(find(ismember(Bands_to_sharpen,Band_overlap_PAN)));
if length(Bands_to_sharpen)~=2^(nextpow2(length(Bands_to_sharpen)))
    disp('Q2n index is advised to work with 2^n bands (where n is an integer)');
end

%% Loading MTF gains at Nyquist frequencies
GNyq_PAN=load_MTF_PAN(sensor_PAN);
GNyq_MS=load_MTF(sensor_MS,im_tag,Bands_to_sharpen_MS);
GNyq=load_MTF(sensor,im_tag,Bands_to_sharpen);

%% PAN and MS scale ratio determination
ratio1_PAN=ratio*ratio_PAN/ratio_MS;
ratio1_MS=ratio;

%% Ground Truth and MS image
if ~exist('ratio_tgt_MS','var'), ratio_tgt_MS=ratio_MS; end
ratio_actual_MS=(size(I_PAN_loaded,1)-2*edge_cut_MS)/(size(I_MS_loaded,1)-2*edge_cut);
flag_imresize_RR=flag_imresize_HS_RR;
I_GT=I_MS_loaded(edge_cut+1:end-edge_cut,edge_cut+1:end-edge_cut,:);
[I_MS_loaded,I_PAN_loaded,edge_temp]=size_check(...
[ratio_tgt_REF,ratio_tgt_REF*ratio_actual_MS/ratio1_MS],...
I_MS_loaded,I_PAN_loaded,'edge',[edge_cut,edge_cut_MS]);
edge_cut=edge_temp(1); edge_cut_MS=edge_temp(2);
[I_MS_LR_e,edge_cut]=scale_image(I_MS_loaded,1/ratio_tgt_REF,flag_interpolation,...
    flag_imresize_RR,'GNyq',GNyq,'edge',edge_cut,'removeedge',0);

%% PAN and MS resizing
ratio1_tgt_MS=min(ratio_actual_MS,ratio_tgt_MS);
[I_PAN,edge_cut_MS]=scale_image(I_PAN_loaded,ratio1_tgt_MS/ratio_actual_MS,...
    flag_interpolation,flag_imresize_MS_sim,'GNyq',GNyq_MS,...
    'edge',edge_cut_MS,'removeedge',0);
[I_PAN,edge_cut_MS]=scale_image(I_PAN,...
    (ratio_tgt_MS/ratio1_tgt_MS)/ratio_tgt_REF,flag_interpolation,...
    flag_imresize_MS_RR,'GNyq',GNyq_MS,'edge',...
    edge_cut_MS,'removeedge',0);
if flag_avoidPANresize~=1
    I_PAN=scale_image(I_PAN,ratio1_MS/ratio_tgt_MS,flag_interpolation,...
        flag_imresize_MS,'GNyq',GNyq_MS,'edge',edge_cut_MS);
    if ratio1_tgt_MS~=ratio1_MS, flag_resizedoriginal_MS=1; else flag_resizedoriginal_MS=0; end
else
    I_PAN=I_PAN(edge_cut_MS+1:end-edge_cut_MS,edge_cut_MS+1:end-edge_cut_MS,:);
    flag_resizedoriginal_MS=0;
end

%% Upsampling
time_EXP=0;
if  flag_avoidEXP~=1
    [I_MS,~,time_EXP]=scale_image(I_MS_LR_e,ratio,flag_interpolation,...
        [],'edge',edge_cut);
end
I_MS_LR=I_MS_LR_e(edge_cut+1:end-edge_cut,edge_cut+1:end-edge_cut,:);

clear('I_HS_LR','I_PAN_temp','I_MS_loaded','I_PAN_loaded',...
    'edge_cut_MS','hshift_PAN','vshift_PAN','hshift_MS',...
    'vshift_MS','hcut','vcut','PAN_choice','MS_choice',...
    'hcut_FR','vcut_FR','Bands_to_sharpen_HS','t2','ii','im_prepare',...
    'flag_1T','flag_HS','start','H','Band_overlap_temp',...
    'sigma','Lfilter','I_MS_temp','filename','edge_temp',...
    'flag_loaddefaultratio','ratio1_PAN','ratio1_MS','ratio_actual_PAN',...
    'ratio_actual_MS');
        
disp('Done');