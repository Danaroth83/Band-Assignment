%% FUSION ALGORITHMS (HS/MS Fusion)

if flag_grouping==7 || strcmpi(flag_grouping,'HYSH'), Fusion_Algorithms_Hypersharpening; return; end

if ~exist('methods_list','var'), methods_list=load_methodslist('custom'); end
if ~exist('flag_thvalues','var'), flag_thvalues=0; end
if ~exist('Bands_to_sharpen','var'), Bands_to_sharpen=1:size(I_MS_LR,3); end
if ~exist('I_GT','var'), flag_FR=1; else flag_FR=0; end
if ~exist('flag_degradation','var'), flag_degradation=0; end
if ~exist('flag_equalization','var'), flag_equalization=0; end
if ~exist('ratio_tgt_REF','var'), ratio_tgt_REF=ratio; end
if ~exist('flag_interpolation','var'), flag_interpolation=0; end
if ~exist('flag_imresize_RR','var'), flag_imresize_RR=0; end
if ~exist('GNyq_MS','var'), GNyq_MS=load_MTF(sensor_PAN,im_tag,Bands_to_sharpen_MS); end
if ~exist('GNyq','var'), GNyq=load_MTF(sensor, im_tag, Bands_to_sharpen); end
if ~exist('flag_grouping','var'), flag_grouping=0; end
if ~exist('shortname','var'), shortname='M'; end

%% Assignation protocol
if strcmpi(shortname,'B') && ismember(flag_grouping,[3:4,20:23])
    error('Unsupported grouping method for this fusion protocol');
end
[Band_grouped,time_grouping]=Band_assigner(flag_grouping,I_PAN,I_MS_LR,ratio,I_MS,GNyq_MS,Band_overlap_MS,Bands_to_sharpen,Bands_to_sharpen_MS,sensor,sensor_PAN,im_tag);
fprintf('Grouping time: %.2f [sec]\n',time_grouping);

%% Fix for unsharpened HS bands (only if flag_grouping==0)
idx_bts=find(ismember(1:length(Bands_to_sharpen),sort(unique(cell2mat(Band_grouped)))));
Bands_to_sharpen=Bands_to_sharpen(idx_bts);
GNyq=GNyq(idx_bts);
I_MS_LR=I_MS_LR(:,:,idx_bts);
I_MS=I_MS(:,:,idx_bts);
if flag_FR==0, I_GT=I_GT(:,:,idx_bts); end
[~,Bands_to_display]=min(abs(repmat((find(idx_bts))',[1,3])-repmat(Bands_to_display,[length(Bands_to_sharpen),1])),[],1);
Band_overlap_MS=cellfun(@(x) intersect(x,idx_bts),Band_overlap_MS,'Un',0);
Band_overlap_MS=cellfun(@(x) find(ismember(idx_bts,x)),Band_overlap_MS,'Un',0);
Band_grouped=cellfun(@(x) changem1(x,1:length(idx_bts),idx_bts), Band_grouped,'Un',0);
% Band_overlap_MS2=sort(unique(cell2mat(Band_overlap_MS)));

%% Removing empty assignments
I_PAN_original=I_PAN;
idx_grouping=~cellfun(@isempty,Band_grouped);
I_PAN=I_PAN(:,:,idx_grouping);
GNyq_MS=GNyq_MS(idx_grouping);

%% Assignation groups
idx_in=Band_grouped(idx_grouping);
Nb_MS=size(I_PAN,3);
Nb_to_sharpen=length(cell2mat(idx_in));
lengths_grp=cumsum(cellfun(@length,idx_in));
idx_out=arrayfun(@(x,y) x:y,[1,lengths_grp(1:end-1)+1],lengths_grp,'Un',0);

%% Low-resolution PAN for equalization
I_PAN_LR=zeros(size(I_PAN));
I_PAN_eq=cell(Nb_MS);
if flag_FR==0
    I_GT_LR=scale_image(I_GT,size(I_PAN,1)/size(I_GT,1),flag_interpolation,...
            flag_imresize_RR,'GNyq',GNyq);
end
for ii=1:Nb_MS
    if flag_FR==0
        [ I_PAN_LR(:,:,ii), I_PAN_eq{ii} ] = PAN_equalization( I_PAN(:,:,ii), I_MS(:,:,idx_in{ii}), flag_degradation, flag_equalization, ratio, GNyq(idx_in{ii}), GNyq_MS(ii),I_GT_LR(:,:,idx_in{ii}));
    else
        [ I_PAN_LR(:,:,ii), I_PAN_eq{ii} ] = PAN_equalization( I_PAN(:,:,ii), I_MS(:,:,idx_in{ii}), flag_degradation, flag_equalization, ratio, GNyq(idx_in{ii}), GNyq_MS(ii));
    end
end

%% Data Fusion
[L1,L2,Nb]=size(I_MS);
NumOutputs=numel(methods_list);
MatrixImage=zeros(L1,L2,Nb_to_sharpen,NumOutputs);
time=zeros(NumOutputs,1);
methods_list_idx=1:NumOutputs;
titleImages=strrep(methods_list,'_','-');
cd Fusion
for ii=1:NumOutputs
    method=methods_list{ii};
    switch method
        case 'EXP'
            for i1=1:Nb_MS
                MatrixImage(:,:,idx_out{i1},ii)=I_MS(:,:,idx_in{i1});
            end
            time(ii) = time_EXP;
        case 'PCA'
            cd CS_new
            t2=tic;
            for i1=1:Nb_MS
                MatrixImage(:,:,idx_out{i1},ii) = PCA_mod(I_MS(:,:,idx_in{i1}),I_PAN(:,:,i1),I_PAN_LR(:,:,i1));
            end
            time(ii) = toc(t2);
            cd ..
        case 'IHS'
            cd CS_new
            t2=tic;
            for i1=1:Nb_MS
                MatrixImage(:,:,idx_out{i1},ii) = IHS_mod(I_MS(:,:,idx_in{i1}),I_PAN(:,:,i1),I_PAN_LR(:,:,i1));
            end
            time(ii) = toc(t2);
            cd ..         
        case 'Brovey'
            cd CS_new
            t2=tic;
            for i1=1:Nb_MS
                MatrixImage(:,:,idx_out{i1},ii) = Brovey_mod(I_MS(:,:,idx_in{i1}),I_PAN(:,:,i1),I_PAN_LR(:,:,i1));
            end
            time(ii) = toc(t2);
            cd ..
        case 'BDSD'
            gcd_temp=gcd(L1,L2);
            if gcd_temp>=128
                BDSD_blocksize=127+find(rem(gcd_temp,128:gcd_temp)==0,1);
            elseif gcd_temp>=32
                BDSD_blocksize=gcd_temp;
            else
                BDSD_blocksize=0;
            end
            if BDSD_blocksize>=Nb+1
                cd BDSD_new
                t2=tic;
                for i1=1:Nb_MS  
                    MatrixImage(:,:,idx_out{i1},ii) = BDSD(I_MS(:,:,idx_in{i1}),I_PAN(:,:,i1),ratio,BDSD_blocksize,GNyq_MS(i1),GNyq(idx_in{i1}));
                end
                time(ii) = toc(t2);
                cd ..
            else
                methods_list_idx=methods_list_idx(methods_list_idx~=ii);
                fprintf('BDSD was suppressed (try using G-BDSD instead)\n');
            end
        case 'GS'
            cd CS_new
            t2=tic;
            for i1=1:Nb_MS
                MatrixImage(:,:,idx_out{i1},ii) = GS_mod(I_MS(:,:,idx_in{i1}),I_PAN(:,:,i1),I_PAN_LR(:,:,i1));
            end
            time(ii) = toc(t2);
            cd ..
        case 'GSA'
            cd GS
            t2=tic;
            for i1=1:Nb_MS
                MatrixImage(:,:,idx_out{i1},ii) = GSA(I_MS(:,:,idx_in{i1}),I_PAN(:,:,i1),I_MS_LR(:,:,idx_in{i1}),ratio);
            end
            time(ii) = toc(t2);
            cd ..
        case 'GSA_local'
            cd GS_new
            t2=tic;
            for i1=1:Nb_MS
                MatrixImage(:,:,idx_out{i1},ii) = GSA_local(I_MS(:,:,idx_in{i1}),I_PAN(:,:,i1),I_MS_LR(:,:,idx_in{i1}),ratio);
            end
            time(ii) = toc(t2);
            cd ..
        case 'PRACS'
            cd PRACS
            t2=tic;
            for i1=1:Nb_MS
                MatrixImage(:,:,idx_out{i1},ii) = PRACS(I_MS(:,:,idx_in{i1}),I_PAN(:,:,i1),ratio);
            end
            time(ii) = toc(t2);
            cd ..
        case 'HPF'
            cd HPF
            t2=tic;
            for i1=1:Nb_MS
                MatrixImage(:,:,idx_out{i1},ii) = HPF(I_MS(:,:,idx_in{i1}),I_PAN(:,:,i1),ratio);
            end
            time(ii) = toc(t2);
            cd ..
        case 'SFIM'
            cd CS_new
            t2=tic;
            for i1=1:Nb_MS
                MatrixImage(:,:,idx_out{i1},ii) = SFIM_mod(I_MS(:,:,idx_in{i1}),I_PAN(:,:,i1),ratio,I_PAN_eq{i1});
            end
            time(ii) = toc(t2);
            cd ..
        case 'Indusion'
            cd Indusion
            t2=tic;
            for i1=1:Nb_MS
                MatrixImage(:,:,idx_out{i1},ii) = Indusion(I_PAN(:,:,i1),I_MS_LR(:,:,idx_in{i1}),ratio);
            end
            time(ii) = toc(t2);
            cd ..
        case 'ATWT'
            cd Wavelet
            t2=tic;
            for i1=1:Nb_MS
                MatrixImage(:,:,idx_out{i1},ii) = ATWT(I_MS(:,:,idx_in{i1}),I_PAN(:,:,i1),ratio);
            end
            time(ii) = toc(t2);
            cd ..
        case 'ATWT_mod'
            cd Wavelet_new
            t2=tic;
            for i1=1:Nb_MS
                MatrixImage(:,:,idx_out{i1},ii) = ATWT_mod(I_MS(:,:,idx_in{i1}),I_PAN(:,:,i1),ratio,I_PAN_eq{i1});
            end
            time(ii) = toc(t2);
            cd ..
        case 'AWLP'
            cd Wavelet
            t2=tic;
            for i1=1:Nb_MS
                MatrixImage(:,:,idx_out{i1},ii) = AWLP(I_MS(:,:,idx_in{i1}),I_PAN(:,:,i1),ratio);
            end
            time(ii) = toc(t2);
            cd ..
        case 'ATWT_M2'
            cd Wavelet
            t2=tic;
            for i1=1:Nb_MS
                MatrixImage(:,:,idx_out{i1},ii) = ATWT_M2(I_MS(:,:,idx_in{i1}),I_PAN(:,:,i1),ratio);
            end
            time(ii) = toc(t2);
            cd ..
        case 'ATWT_M3'
            cd Wavelet
            t2=tic;
            for i1=1:Nb_MS
                % MatrixImage(:,:,idx_out{i1},ii) = ATWT_M3(I_MS(:,:,idx_in{i1}),I_PAN(:,:,i1),ratio); 
                %%% There was no warning for complex numbers in RR, but I changed to real regardless
                MatrixImage(:,:,idx_out{i1},ii) = real(ATWT_M3(I_MS(:,:,idx_in{i1}),I_PAN(:,:,i1),ratio));
            end
            time(ii) = toc(t2);
            cd ..
        case 'MTF_GLP'
            cd GLP_new
            t2=tic;
            for i1=1:Nb_MS
                MatrixImage(:,:,idx_out{i1},ii) = MTF_GLP(I_PAN(:,:,i1),I_MS(:,:,idx_in{i1}),GNyq(idx_in{i1}),ratio);
                % MatrixImage(:,:,idx_out{i1},ii) = MTF_GLP_ECB(I_MS(:,:,idx_in{i1}),I_PAN(:,:,i1),ratio,[9 9],2.5,GNyq(idx_in{i1}));
                % MatrixImage(:,:,idx_out{i1},ii) = MTF_GLP_CBD(I_MS(:,:,idx_in{i1}),I_PAN(:,:,i1),ratio,[55 55],-Inf,GNyq(idx_in{i1}));
            end
            time(ii) = toc(t2);
            cd ..
        case 'MTF_GLP_HPM_PP'
            if 2^nextpow2(ratio)==ratio,
                cd GLP_new
                t2=tic;
                for i1=1:Nb_MS
                    MatrixImage(:,:,idx_out{i1},ii) = MTF_GLP_HPM_PP(I_PAN(:,:,i1),I_MS_LR(:,:,idx_in{i1}),GNyq(idx_in{i1}),ratio);
                end
                time(ii) = toc(t2);
                cd ..
            else
                methods_list_idx=methods_list_idx(methods_list_idx~=ii);
                fprintf('MTF-GLP-HPM With PostProcessing was suppressed\n');
            end
        case 'MF_HG'
            individuals=[1 0 1 0 0 0]; %unused
            % ReDefines also operator after MF_Adapt_Settings 1st and 2nd bits:
            operators = [0 0 0];
            % 000: gradient; 010: top-hat; 100: toggle; 110: Bai (toggle+top-hat)
            % 001: LOCO; 011:median; 101: reconstuction
            textse=repmat([ 0 1 0 1 1 1 0 1 0 0 1],1,3);
            cd MF_Adapt
            save chromosomes_file textse individuals operators
            t2=tic;
            for i1=1:Nb_MS
                if flag_FR==0
                    [MatrixImage(:,:,idx_out{i1},ii),dum,D_MF_HG,P_LP_MF_HG]=MF_Adapt_gradLoad(I_PAN(:,:,i1),I_MS(:,:,idx_in{i1}),L,ratio,sensor,I_GT_LR(:,:,idx_in{i1}));
                else
                    [MatrixImage(:,:,idx_out{i1},ii),dum,D_MF_HG,P_LP_MF_HG]=MF_Adapt_gradLoad(I_PAN(:,:,i1),I_MS(:,:,idx_in{i1}),L,ratio,sensor,[]);%whole
                end
            end
            time(ii) = toc(t2);
            cd ..
        case 'MTF_GLP_HPM'
            cd GLP_new
            t2=tic;
            for i1=1:Nb_MS
                % MatrixImage(:,:,idx_out{i1},ii) = MTF_GLP_HPM(I_PAN(:,:,i1),I_MS(:,:,idx_in{i1}),GNyq(idx_in{i1}),ratio);
                % MatrixImage(:,:,idx_out{i1},ii) = MTF_GLP_HPM2(I_PAN(:,:,i1),I_MS(:,:,idx_in{i1}),GNyq(idx_in{i1}),ratio);  %Tryouts
                MatrixImage(:,:,idx_out{i1},ii) = MTF_GLP_HPM_mod(I_PAN(:,:,i1),I_MS(:,:,idx_in{i1}),GNyq(idx_in{i1}),ratio,I_PAN_eq{i1});
            end
            time(ii) = toc(t2);
            cd ..
        case 'MTF_GLP_CBD'
            cd GS
            t2=tic;
            for i1=1:Nb_MS
                MatrixImage(:,:,idx_out{i1},ii) = GS2_GLP_mod(I_MS(:,:,idx_in{i1}),I_PAN(:,:,i1),ratio,GNyq(idx_in{i1}));
            end
            time(ii) = toc(t2);
            cd ..
        case 'GFPCA'
            No_PCs=4;    % number of principal components used by guided filtering with pan/RGB image, e.g. No_PCs=5
            r=8;         % the size of local sliding window, e.g. r=8
            eps=0.001^2; % regularization parameter determining the degree of blurring for the guided filter, e.g. eps = 0.01^2        
            cd GFPCA
            t2=tic;
            for i1=1:Nb_MS
                MatrixImage(:,:,idx_out{i1},ii) = GFPCA(I_MS_LR(:,:,idx_in{i1}), I_PAN(:,:,i1), min(No_PCs,length(idx_in{i1})), r, eps);
            end
            time(ii) = toc(t2);
            cd ..
        case 'CNMF'
            cd CNMF
            t2=tic;
            % for i1=1:Nb_MS
            %     MatrixImage(:,:,idx_out{i1},ii) = CNMF_fusion(I_MS_LR(:,:,idx_in{i1}),I_PAN(:,:,i1));
            % end
             MatrixImage_temp = CNMF_fusion(I_MS_LR,I_PAN_original);
            for i1=1:Nb_MS
                MatrixImage(:,:,idx_out{i1},ii) = MatrixImage_temp(:,:,idx_in{i1});
            end
            time(ii) = toc(t2);
            cd ..
        case {'BayesNaive','BayesSparse'}
            if strcmpi(method,'BayesNaive'), prior='Gaussian'; end
            if strcmpi(method,'BayesSparse'), prior='Sparse'; end
            cd ..
            [KerBlu,start_pos]=load_KerBlu(im_tag,ratio,GNyq(idx_in{i1}),min(size(I_PAN,1),size(I_PAN,2)));
            cd 'Fusion\BayesFusion'
            t2=tic;
            setup;           
            MatrixImage_temp=BayesianFusion_mod2(I_MS_LR,I_PAN_original,Band_overlap_MS,KerBlu,ratio,prior,start_pos);
            for i1=1:Nb_MS
                MatrixImage(:,:,idx_out{i1},ii) = MatrixImage_temp(:,:,idx_in{i1});
            end
            time(ii) = toc(t2);
            cd ..
        case {'BayesNaive_new','BayesSparse_new','BayesML','BayesNLM'}
            if strcmpi(method,'BayesNaive_new'), prior='Gaussian'; end
            if strcmpi(method,'BayesSparse_new'), prior='Sparse'; end
            if strcmpi(method,'BayesML'), prior='ML'; end
            if strcmpi(method,'BayesNLM'), prior='NLM'; end
            Band_overlap_MS_temp=Band_overlap_MS;
            for i1=1:numel(Band_overlap_MS)
                if length(Band_overlap_MS_temp{i1})<=1
                    a1=sort(unique(cell2mat(Band_overlap_MS(1:max(1,i1-1)))));
                    a2=sort(unique(cell2mat(Band_overlap_MS(min(numel(Band_overlap_MS),i1+1):end))));
                    if isempty(a1)
                        Band_overlap_MS_temp{i1}=a2(1:2);
                    elseif isempty(a2)
                        Band_overlap_MS_temp{i1}=a1(end-1:end);
                    else
                        Band_overlap_MS_temp{i1}=[a1(end),a2(1)];
                    end
                end
            end
            cd ..
            [KerBlu,start_pos]=load_KerBlu(im_tag,ratio,GNyq,min(size(I_PAN_original,1),size(I_PAN_original,2)));
            cd 'Fusion\BlindFuse-master'
            t2=tic;           
            MatrixImage_temp=main_BlindFuse_mod(I_PAN_original,I_MS_LR,ratio,Band_overlap_MS_temp,start_pos,prior,'PCA');
            for i1=1:Nb_MS
                MatrixImage(:,:,idx_out{i1},ii) = MatrixImage_temp(:,:,idx_in{i1});
            end
            time(ii) = toc(t2);
            cd ..
        case 'HySure'
            cd ..
            [~,start_pos]=load_KerBlu(im_tag,ratio,GNyq(idx_in{i1}),min(size(I_PAN,1),size(I_PAN,2)));
            cd 'Fusion\HySure'
            t2=tic;
            % for i1=1:Nb_MS
            %     MatrixImage(:,:,idx_out{i1},ii)= HySure_wrapper_mod(I_MS_LR(:,:,idx_in{i1}), I_PAN(:,:,i1), ratio, 1:length(idx_in{i1}));
            % end
            %% Temporary fix for empty overlapping bands
            Band_overlap_MS_temp=Band_overlap_MS;
            for i1=1:numel(Band_overlap_MS)
                if length(Band_overlap_MS_temp{i1})<=1
                    a1=sort(unique(cell2mat(Band_overlap_MS(1:max(1,i1-1)))));
                    a2=sort(unique(cell2mat(Band_overlap_MS(min(numel(Band_overlap_MS),i1+1):end))));
                    if isempty(a1)
                        Band_overlap_MS_temp{i1}=a2(1:2);
                    elseif isempty(a2)
                        Band_overlap_MS_temp{i1}=a1(end-1:end);
                    else
                        Band_overlap_MS_temp{i1}=[a1(end),a2(1)];
                    end
                end
            end
            MatrixImage_temp=HySure_wrapper_MS(I_MS_LR,I_PAN_original,ratio,Band_overlap_MS_temp,start_pos,im_tag);
            for i1=1:Nb_MS
                MatrixImage(:,:,idx_out{i1},ii) = MatrixImage_temp(:,:,idx_in{i1});
            end
            time(ii) = toc(t2);
            cd ..
        case 'G_BDSD'
            cd BDSD_new
            t2=tic;
            for i1=1:Nb_MS
                MatrixImage(:,:,idx_out{i1},ii) = GBDSD_Fusion(I_PAN(:,:,i1),I_MS(:,:,idx_in{i1}),ratio,GNyq_MS(i1),GNyq(idx_in{i1}));
            end
            time(ii) = toc(t2);
            cd ..
        case 'L_BDSD'
            cd BDSD_new
            t2=tic;
            for i1=1:Nb_MS
                MatrixImage(:,:,idx_out{i1},ii) = IBDSD_Fusion(I_PAN(:,:,i1),I_MS(:,:,idx_in{i1}),ratio,GNyq_MS(i1),GNyq(idx_in{i1}));
                % MatrixImage(:,:,idx_out{i1},ii) = LBDSD_Fusion(I_PAN(:,:,i1),I_MS(:,:,idx_in{i1}),ratio,GNyq_MS(i1),GNyq(idx_in{i1}));
            end
            time(ii) = toc(t2);
            cd ..
        case 'C_BDSD'
            cd BDSD_new
            t2=tic;
            [MatrixImage(:,:,idx_out{1},ii), index_map_FR, index_map, Nclusters]=CBDSD_Fusion(I_PAN(:,:,1),I_MS(:,:,idx_in{1}),ratio,GNyq_MS(1),GNyq(idx_in{1}));
            for i1=2:Nb_MS
                MatrixImage(:,:,idx_out{i1},ii) = CBDSD_Fusion(I_PAN(:,:,i1),I_MS(:,:,idx_in{i1}),ratio,GNyq_MS(i1),GNyq(idx_in{i1}),index_map_FR,index_map,Nclusters);
            end
            % [MatrixImage(:,:,idx_out{1},ii), index_map_FR, Nclusters] = SCBDSD_Fusion(I_PAN(:,:,1),I_MS(:,:,idx_in{1}),ratio,GNyq_MS(1),GNyq(idx_in{1}));
            % for i1=2:Nb_MS
            %     MatrixImage(:,:,idx_out{i1},ii) = SCBDSD_Fusion(I_PAN(:,:,i1),I_MS(:,:,idx_in{i1}),ratio,GNyq_MS(i1),GNyq(idx_in{i1}),index_map_FR,Nclusters);
            % end
            time(ii) = toc(t2);
            cd ..
       case 'DMSC'
            cd Demosaicking

            cd ../GS
            spectralweights=estimation_alpha(I_MS_LR,imresize(I_PAN,1/ratio),'global');
            cd ../Demosaicking
            I_PAN_comp=zeros(size(I_PAN));
            t2=tic;
            for i1=1:Nb_MS
                [MatrixImage(:,:,idx_out{i1},ii),I_PAN_comp(:,:,i1)] = main4(I_PAN(:,:,i1),I_MS_LR(:,:,idx_in{i1}),ratio,GNyq(idx_in{i1}),spectralweights(idx_in{i1}),lambda_DMSC);
            end
            time(ii) = toc(t2);
            cd ..
    end
    if ismember(ii,methods_list_idx)
        fprintf('Elaboration time %s: %.2f [sec]\n',titleImages{ii},time(ii));
    end
end
cd ..

%% Removal of failed fusions
methods_list=methods_list(methods_list_idx);
titleImages=titleImages(methods_list_idx);
MatrixImage=MatrixImage(:,:,:,methods_list_idx);
time=time(methods_list_idx);
NumOutputs=length(methods_list_idx);

%% Reordering image matrix
MatrixImage_in=MatrixImage;
MatrixImage=[];
idx_in=cell2mat(idx_in);
for ii=1:Nb
    idx=ismember(idx_in,ii);
    MatrixImage=cat(3,MatrixImage,mean(MatrixImage_in(:,:,idx,:),3));
    MatrixImage_in=MatrixImage_in(:,:,~idx,:);
    idx_in=idx_in(~idx);
end

%% Upsampling of fused images
MatrixImage_in=MatrixImage;
MatrixImage=[];
for ii=1:size(MatrixImage_in,4)
    MatrixImage_temp=scale_image(MatrixImage_in(:,:,:,1),ratio_tgt_REF/ratio,...
    flag_interpolation,flag_imresize_RR,'GNyq',GNyq);
    MatrixImage_in(:,:,:,1)=[];
    MatrixImage=cat(4,MatrixImage,MatrixImage_temp);
end

%% Cuts values outside of radiometric range
if flag_thvalues==1
    MatrixImage(MatrixImage<2^L)=2^L;
    MatrixImage(MatrixImage<0)=0;
end

%% Restoring original values
I_PAN=I_PAN_original;

clear('I_PAN_LR','I_PAN_eq','I_GT_LR','ii','L1','L2','Nb','NumOutputs',...
    'methods_list_idx','t2','time_EXP','GNyq','GNyq_MS',...
    'BDSD_blocksize','dum','D_MF_HG','P_LP_MF_HG','individuals',...
    'operators','textse','No_PCs','r','eps','KerBlu','start_pos',...
    'MatrixImage_temp','I_MS_LR_temp','Band_overlap_temp',...
    'index_map_FR', 'index_map', 'Nclusters','MatrixImage_in','idx_in',...
    'idx_out','idx_bts','a1','a2','lengths_grp','I_GT_LR')