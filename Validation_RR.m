%% REDUCED RESOLUTION VALIDATION

%% Print Options
if ~exist('printEPS','var'), printEPS = 0; end
if ~exist('printTEX','var'), printTEX = 1; end
if ~exist('printTIFF','var'), printTIFF = 0; end
if ~exist('printMAT','var'), printMAT = 1; end
if ~exist('flag_draw','var'), flag_draw = 1; end

%% Cut Final Image
if ~exist('flag_cutbounds_vis','var'), flag_cutbounds_vis = 0; end % Removes edges in visualization
if ~exist('flag_cutbounds_val','var'), flag_cutbounds_val = 1; end % Removes edges in validation
if ~exist('dim_cut','var'), dim_cut = 11; end                      % Amount of edge removal in visualization
if ~exist('flag_thvalues_vis','var'), flag_thvalues_vis = 0; end   % Threshold values out of dynamic range during visualization
if ~exist('flag_thvalues_val','var'), flag_thvalues_val = 0; end   % Threshold values out of dynamic range during validation

if ~exist('sensor_HS_list','var'), sensor_HS_list={'HYP','AVIRIS','CHR','ROSIS'}; end
if ~exist('flag_resizedoriginal_PAN','var'), flag_resizedoriginal_PAN=0; end
if ~exist('flag_resizedoriginal_MS','var'), flag_resizedoriginal_MS=0; end
if ~exist('flag_RR_val','var'), flag_RR_val=1; end
if ~exist('flag_R2','var'), flag_R2 = 0; end
if ~exist('shortname','var')
    if size(I_PAN,3)==1, shortname='P'; else shortname='M'; end
end

%% Fix for full resolution validation
if flag_RR_val==5
    flag_RRtoFR=1;
    Validation_FR;
    return;
end

if any(strcmp(sensor,sensor_HS_list)), lowres_string='HS'; else lowres_string='MS'; end
if any(strcmp(shortname,{'P','S'})), hires_string='PAN'; else hires_string='MS'; end
if flag_noautoblock==1, Qblocks_size=Qblocks_size_setup; end

%% Image visualization and printing
if strcmpi(shortname,{'P'})
    filename_output=sprintf('%s_%s_p%df%dr%d',im_tag,shortname,round(ratio_tgt_PAN),round(ratio),round(ratio_tgt_REF));
elseif strcmpi(shortname,{'M'})
    filename_output=sprintf('%s_%s_m%df%dr%d',im_tag,shortname,round(ratio_tgt_MS),round(ratio),round(ratio_tgt_REF));
else
    filename_output=sprintf('%s_%s_p%dm%df%dr%d',im_tag,shortname,round(ratio_tgt_PAN),round(ratio_tgt_MS),round(ratio),round(ratio_tgt_REF));
end
if any(strcmpi(shortname,{'P','S'}))
    filename_output=[filename_output,'_RR'];
else
   flag_grouping_list=[0,1,2,3,4,5,6,7,8,9,10,12,14,15,20,21,22,23,24,25,26,27];
   grouping_label_list={'OVLP','SAM','CC','MSD','SAM2','SAMHS','SAMHSEXP','HYSH','CCEXP','SAMHSEXP2','CCHM','CCFR','SAMFR','SAMHSFR','CEN','DIS','CEN2','DIS2','CCOV','RSR','LSQ','LSQ2'};
   grouping_label=grouping_label_list{flag_grouping_list==flag_grouping};
   filename_output=sprintf([filename_output,'_g%s_RR'],grouping_label);
end
[L1,L2,Nb,NumOutputs]=size(MatrixImage);

if flag_draw~=0
    Xmap=10; Ymap=19;
    text_filename=[{[lowres_string,'_LR'],hires_string,'GT'},methods_list];
    text_caption=[{[lowres_string,'-LR'],hires_string,'GT'},titleImages];
    cd Visualization
    for ii=-2:NumOutputs
        idx=ii+3;
        if idx==1
            showImageHSLR(I_MS_LR,printEPS,idx,flag_cutbounds_vis,dim_cut,flag_thvalues_vis,L,Bands_to_display,ratio);
        elseif idx==2
            if any(strcmpi(shortname,{'P','S'}))
                showPan(I_PAN,printEPS,idx,flag_cutbounds_vis,dim_cut);
            else
                showImageHS(I_PAN,printEPS,idx,flag_cutbounds_vis,dim_cut,0,[],Bands_to_display_MS);
            end
        elseif idx==3
            showImageHS(I_GT,printEPS,idx,flag_cutbounds_vis,dim_cut,flag_thvalues_vis,L,Bands_to_display);
        else
            showImageHS(MatrixImage(:,:,:,ii),printEPS,idx,flag_cutbounds_vis,dim_cut,flag_thvalues_val,L,Bands_to_display);
        end
        if printTIFF==1
            filename_TIFF=['../Outputs/',filename_output,'_',text_filename{idx},'.tif'];
            option_plot=get(gcf,'PaperPositionMode');
            set(gcf,'PaperPositionMode','auto');
            print(filename_TIFF,'-dtiffn','-r0');
            set(gcf,'PaperPositionMode',option_plot);
        end
        if printEPS==1
            close(idx+1);
            filename_EPS=['../Outputs/',filename_output,'_',text_filename{idx},'.eps'];
            movefile(['../Outputs/',num2str(idx),'.eps'],filename_EPS);
        end
        text(Xmap,Ymap,text_caption{idx},'EdgeColor','k','BackgroundColor','w');
    end
    drawnow;
    cd ..
else
    if printTIFF==1, disp('TIFF images weren''t saved'); end
    if printEPS==1, disp('EPS images weren''t saved'); end
end

%% Run validation on a subset of bands
if any(strcmpi(shortname,{'P','S'}))
    Band_overlap=Band_overlap_PAN;
else
    Band_overlap=sort(unique(cell2mat(Band_overlap_MS)));
end
bandselect_val=1:Nb;
if flag_bandselect_val==1
    bandselect_val=Band_overlap;
elseif flag_bandselect_val==2
    bandselect_val=bandselect_val(~ismember(bandselect_val,Band_overlap));
end
if isempty(bandselect_val), error('No band was selected for validation'); end
if flag_bandselect_val~=0
    save('Output_temp.mat','I_MS_LR','I_GT','MatrixImage');
    I_MS_LR=I_MS_LR(:,:,bandselect_val);
    I_GT=I_GT(:,:,bandselect_val);
    MatrixImage=MatrixImage(:,:,bandselect_val,:);
end

%% Evaluation of quality indexes

if flag_bandbyband==1
    SCC_b=zeros(Nb,NumOutputs);
    Q_b=zeros(Nb,NumOutputs);
    ERGAS_b=zeros(Nb,NumOutputs);
end
if flag_RR_val==4
    qindex_list={'SCCo_laplacian','SCCo_sobel','SCCo_sobel2',...
        'SCCo_prewitt','SCCo_prewitt2','SCC_laplacian',...
        'SCC_sobel','SCC_sobel2','SCC_prewitt','SCC_prewitt2'};
elseif flag_RR_val==3
    qindex_list={'Q2n_b16s16','Q2^n_b32s32','Q2n_b32s8','Q2n_b36s36',...
        'Q2n_b32s32c1','Q2n_b32s32c10','Q2n_alt_b36s36','Q2n_alt_b36s36'};
elseif flag_RR_val==2
    qindex_list={'Q2^n','SAM','ERGAS','SCC','UIQI','Q2n_alt'};
elseif flag_RR_val==1
    qindex_list={'Q2^n','SAM','ERGAS','SCC','UIQI'};
end
NumIndexes=numel(qindex_list);
MatrixResults=zeros(NumOutputs,NumIndexes);

fprintf('Calculating quality indices:  0/%2d',NumOutputs);
for ii=1:NumOutputs
    method=methods_list{ii};
    [MatrixResults(ii,:),columnLabels] = indexes_evaluation_mod2(MatrixImage(:,:,:,ii),I_GT,qindex_list,ratio_tgt_REF,L,Qblocks_size,flag_cutbounds_val,dim_cut,flag_thvalues_val);
    if flag_bandbyband==1
        [SCC_b(:,ii),Q_b(:,ii),ERGAS_b(:,ii)]=indexes_evaluation_bandbyband(MatrixImage(:,:,:,ii),I_GT,ratio_tgt_REF);
    end
    fprintf([repmat('\b',[1,5]),'%2d/%2d'],ii,NumOutputs);
end
fprintf('. Done!\n');

%% Print in LaTEX
filename_TEX=['Outputs/',filename_output,sprintf('_b%d_v%dq%d.tex',flag_bandselect_val,flag_RR_val,Qblocks_size)];
switch shortname
    case 'P'
        fullname=sprintf('%s/PAN-%s Fusion',lowres_string,sensor_PAN);
    case 'S'
        fullname=sprintf('%s/(Simulated-PAN)-%s Fusion',lowres_string,sensor_PAN);
    case 'M'
        fullname=sprintf('%s/MS-%s Fusion',lowres_string,sensor_PAN);
    case 'B'
        fullname=sprintf('%s/(PAN+MS)-%s Fusion',lowres_string,sensor_PAN);
    otherwise
        fullname='Fusion results';
end
lowres_resolution=load_resolution(sensor,im_tag,lowres_string)*ratio_tgt_REF;
if printTEX==1
    [L1,L2,Nb]=size(I_GT);
    load('PreambleTEX.mat');
    fid=fopen(filename_TEX,'w');
    fprintf(fid,PreambleTEX);
    fprintf(fid,'\\begin{document}\n\n');
    fprintf(fid,'\\section*{%s\\\\\n%s (Reduced Resolution)}\n',fullname,strrep(im_tag,'_','-'));
    fprintf(fid,'(Bands: %d - Ratio: %d - GT sizes: %d x %d - Sensor: %s+%s)\\\\\n',Nb,ratio_tgt_REF,L2,L1,sensor,sensor_PAN);
    fprintf(fid,'%s-%s sharpening realized at target resolution %.12gm with scale ratio %d\\\\\n',lowres_string,sensor,lowres_resolution/ratio,ratio);
    
    fprintf(fid,'Initial resolution: ');
    if ~strcmpi(shortname,'M')
        fprintf(fid,'PAN-%s - %.12gm',sensor_PAN,lowres_resolution/ratio_tgt_PAN);
    end
    if ~strcmpi(shortname,'P')
        if ~strcmpi(shortname,'M'), fprintf(fid,'; '); end
        fprintf(fid,'MS-%s - %.12gm',sensor_PAN,lowres_resolution/ratio_tgt_MS);
    end
    fprintf(fid,'\\\\\n');
    flag_bandselect_val_list=[0,1,2];
    bandselect_list={'All bands','Overlapping bands','Non overlapping bands'};
    fprintf(fid,'Amount of validated bands: %d (%s)\\\\\n', length(bandselect_val),bandselect_list{flag_bandselect_val_list==flag_bandselect_val});
    if any(strcmpi(shortname,{'M','F'}))
        fprintf(fid,'Amount of MS-%s bands used for fusion: %d\\\\\n',sensor_PAN,size(I_PAN,3));
    elseif strcmpi(shortname,'B')
        fprintf(fid,'Amount of MS-%s bands used for fusion: %d\\\\\n',sensor_PAN,size(I_PAN,3)-1);
    end
    if strcmpi(shortname,'B')
        fprintf(fid,'PAN-%s added to MS-%s as extra band\\\\\n',sensor_PAN,sensor_PAN);
    elseif strcmpi(shortname,'F')
        fprintf(fid,'MS-%s previously sharpened with PAN-%s at resolution %.12gm with scale ratio %d with ATWT\\\\\n',...
            sensor_PAN,sensor_PAN,lowres_resolution/ratio_tgt_PAN,ratio_tgt_PAN/ratio_tgt_MS);
    end
    if any(strcmpi(shortname,{'M','F','B'}))
        flag_grouping_list=[0,1,2,3,4,5,6,7,8,9,10,12,14,15,20,21,22,23,24,25,26,27];
        grouping_label_list={'Band overlap','SAM-based','CC-based','Minimum central wavelength distance',...
            'SAM-based (simplified)', 'SAM-based with simulated HS', 'SAM-based with simulated interpolated HS',...
            'Hypersharpening','sCC-based with interpolated HS','SAM-based with degraded MS kept at its scale',...
            'CC-based with histogram matching','CC-based at FR','SAM based (simplified) at FR','SAM-based with simulated HS at FR',...
            'Estimated centroid band distance','Estimated dispersion-based ML',...
            'Given center band distance','Given dispersion-based ML','CCOV-based',...
            'Radiometric Spatial Ratio-based','Best HS regression weights to simulate MS',...
            'Best MS regression weights to simulate HS'};
        grouping_label=grouping_label_list{flag_grouping_list==flag_grouping};
        fprintf(fid,'%s grouping: %s\\\\\n',lowres_string,grouping_label);
    end
    if flag_R2==1, fprintf(fid,'R2 index: %.4f\\\\\n',R2_idx); end
    flag_imresize_list=[0,1,2];
    downsample_list={'MTFfiltered','imresize','HS-MTFfiltered'};
    if any(strcmpi(shortname,{'P','B','F','S'}))
        fprintf(fid,'PAN: RR downsampling method - %s\\\\\n',downsample_list{flag_imresize_list==flag_imresize_PAN_RR});
        if flag_resizedoriginal_PAN==1
           fprintf(fid,'PAN: Further downsampling method - %s\\\\\n',downsample_list{flag_imresize_list==flag_imresize_PAN});
        end
    end
    if any(strcmpi(shortname,{'M','B','F','S'})) 
        fprintf(fid,'MS: RR downsampling method - %s\\\\\n',downsample_list{flag_imresize_list==flag_imresize_MS_RR});
         if flag_resizedoriginal_MS==1
           fprintf(fid,'MS: Further downsampling method - %s\\\\\n',downsample_list{flag_imresize_list==flag_imresize_MS});
        end
    end
    fprintf(fid,'%s: RR downsampling method - %s\\\\\n',lowres_string,downsample_list{flag_imresize_list==flag_imresize_RR});
    
    flag_interpolation_list=[0,1,2];
    interpol_string={'Hamming-windowed kernel','Bicubic kernel','Lagrange polynomial kernel'};
    fprintf(fid,'Interpolation: %s\\\\\n',interpol_string{flag_interpolation_list==flag_interpolation});
    
    fprintf(fid,'Equalization - PAN Degradation: ');
    flad_degradation_list=[0,1,2,3,4,5];
    degrad_string={'Undegraded','Bicubic','PAN-MTF/matched','ATWT-based','Hamming windowed ideal','MS-MTF/matched'};
    fprintf(fid,degrad_string{flad_degradation_list==flag_degradation});
    fprintf(fid,' - Method: ');
    flag_equalization_list=[0,1,2,3];
    equal_string={'Classic','Classic with GT','Regression based','PCA regression based'};
    fprintf(fid,'%s\\\\\n',equal_string{flag_equalization_list==flag_equalization});
    
    fprintf(fid,'\\\\\n');
    fclose(fid);
    matrix2latex_mod3(MatrixResults,filename_TEX, 'rowLabels',titleImages,...
            'columnLabels',columnLabels,'alignment','c','format', ...
            '%.4f','highlight',matrix2latex_highlightoption(columnLabels),...
            'hCol',ones(1,length(columnLabels)));
    fid=fopen(filename_TEX,'a');
    fprintf(fid,'\n\\end{document}');
    fclose(fid);
    copyfile(filename_TEX,'Output.tex');
end  

if flag_bandselect_val~=0, load('Output_temp.mat'); delete('Output_temp.mat'); end

%% View All

% cd Visualization
% figure, showImagesAll_print(MatrixImage,titleImages,fliplr(Bands_to_display),flag_cutbounds_vis,dim_cut,0,printEPS);
% cd ..
% if printEPS==1
%     for ii=1:NumOutputs, close(NumOutputs+3+ii); end
% end

%% Create easy to read variables and clearing MatrixImage
% for ii=1:NumOutputs
%     assignin('caller',matlab.lang.makeValidName(['time_',methods_list{ii}]),time(ii));
%     assignin('caller',matlab.lang.makeValidName(['I_',methods_list{ii}]),MatrixImage(:,:,:,1));
%     MatrixImage(:,:,:,1)=[];
% end

if printMAT==1
    save([filename_TEX(1:end-4),'.mat'],'MatrixResults','columnLabels',...
        'titleImages','shortname','ratio_tgt_REF','ratio','flag_interpolation',...
        'bandselect_val','fullname','flag_degradation','flag_equalization',...
        'flag_imresize_RR');
    if ~strcmpi(shortname,'M')
        save([filename_TEX(1:end-4),'.mat'],'flag_imresize_PAN',...
            'flag_imresize_PAN_RR','ratio_tgt_PAN','-append');
    end
    if ~strcmpi(shortname,'P')
        save([filename_TEX(1:end-4),'.mat'],'flag_imresize_MS',...
            'flag_imresize_MS_RR','ratio_tgt_MS','-append');
    end
    if ~any(strcmpi(shortname,{'P','S'}))
        save([filename_TEX(1:end-4),'.mat'],'flag_grouping','-append');
    end
    if flag_bandbyband==1
        save([filename_TEX(1:end-4),'.mat'],'SCC_b','ERGAS_b','Q_b','-append');
    end
end

clear('ii','i1','idx','current_flag','method','I_temp','L1','L2','Nb',...
    'Xmap','Ymap','chosen_methods','t2','label','NumOutputs','L_methods_list',...
    'sig','size_kernel','No_PCs','r','eps','flag_all',...
    'filename_TIFF','filename_TEX','PreambleTEX','gcd_temp',...
    'fid','option_plot','text_filename','text_caption',...
    'lowres_string','flag_interpolation_list','interpol_list',...
    'flag_bandselect_val_list','degrad_string','downsample_list',...
    'bandselect_list');

if printMAT==1
    save('Output.mat');
    copyfile('Output.mat',['Outputs/',filename_output,'.mat']);
end