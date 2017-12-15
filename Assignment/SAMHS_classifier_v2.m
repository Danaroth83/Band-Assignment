function [ SAM_groups,SAM_mat ] = SAMHS_classifier_v2( I_PAN,I_MS_LR,ratio,GNyq )
%SAM_CLASSIFIER Classification of HS imagery based on SAM index
%   Grouping of HS bands to fuse with MS imagery;
%   Input images are:
%   I_PAN: MS images
%   I_MS:  HS images
%   ratio: scale ratio between MS and HS images
%   GNyq: Gain at Nyquist frequency
%   band_overlap_cell: cell of bands of HS overlapping with each MS

Nb_MS=size(I_PAN,3);
[L1,L2,Nb]=size(I_MS_LR);

if Nb_MS==1
    SAM_groups={1:Nb};
    return;
end

cd ../Quality_Indices_HS
I_PAN_LR=MTF_filter(I_PAN,GNyq,ratio,41);
I_PAN_LR=I_PAN_LR(floor(ratio/2)+1:ratio:end,floor(ratio/2)+1:ratio:end,:);
mean2_I_PAN_LR=squeeze(sum(sum(I_PAN_LR,1),2))/L1/L2;
std2_I_PAN_LR=sqrt(squeeze(sum(sum(I_PAN_LR.^2,1),2))/L1/L2-mean2_I_PAN_LR.^2);
cd ../Assignment

mean2_I_MS_LR=squeeze(sum(sum(I_MS_LR,1),2))/L1/L2;
std2_I_MS_LR=sqrt(squeeze(sum(sum(I_MS_LR.^2,1),2))/L1/L2-mean2_I_MS_LR.^2);
I_PAN_HM=zeros(L1,L2,Nb_MS,Nb);
for ii1=1:Nb_MS
    for ii2=1:Nb
        I_PAN_HM(:,:,ii1,ii2)=(I_PAN_LR(:,:,ii1)-mean2_I_PAN_LR(ii1)*ones(L1,L2))*std2_I_MS_LR(ii2)/std2_I_PAN_LR(ii1)+mean2_I_MS_LR(ii2)*ones(L1,L2);
    end
end

cd ../Quality_indices_HS
% SAM_mat=zeros(Nb_MS,Nb);
% for ii1=1:Nb_MS
%     for ii2=1:Nb
%         I_MS_simul=I_MS_LR;
%         I_MS_simul(:,:,ii2)=I_PAN_HM(:,:,ii1,ii2);
%         SAM_mat(ii1,ii2)=SAM(I_MS_simul,I_MS_LR);
%     end
% end
SAM_mat=SAM_HSMS(I_MS_LR,I_PAN_HM);
[~,min_index]=min(SAM_mat);

cd ../Assignment
SAM_groups=cell(1,Nb_MS);
for ii1=1:Nb_MS
    SAM_groups{ii1}=find(min_index==ii1);
end
end

