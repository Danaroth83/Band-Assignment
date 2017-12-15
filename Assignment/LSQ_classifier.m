function [ LSQ_groups ] = LSQ_classifier( I_PAN,I_MS )
%CC_CLASSIFIER Grouping HS images for fusing with MS
%   I_MS:    interpolated HS images
%   I_PAN:   original MS image

[L1,L2,Nb]=size(I_MS);
Nb_MS=size(I_PAN,3);

I_MS=reshape(I_MS,[L1*L2,Nb]);
I_PAN=reshape(I_PAN,[L1*L2,Nb_MS]);

alpha=zeros(Nb_MS,Nb);
mean2_I_PAN=sum(I_PAN,1)/L1/L2;
std2_I_PAN=sqrt(sum(I_PAN.^2,1)/L1/L2-mean2_I_PAN.^2);
mean2_I_MS=sum(I_MS,1)/L1/L2;
std2_I_MS=sqrt(sum(I_MS.^2,1)/L1/L2-mean2_I_MS.^2);
for ii=1:Nb
    I_PAN_temp=(I_PAN-repmat(reshape(mean2_I_PAN,[1,Nb_MS]),[L1*L2,1]))./repmat(reshape(std2_I_PAN,[1,Nb_MS]),[L1*L2,1])*std2_I_MS(ii)+repmat(mean2_I_MS(ii),[L1*L2,Nb_MS]);
    alpha(:,ii)=lsqnonneg(I_PAN_temp,I_MS(:,ii));
    % alpha(:,ii)=alpha(:,ii)/sum(alpha(:,ii));
end

[~,groups]=max(alpha);

LSQ_groups=cell(1,Nb_MS);
for ii1=1:Nb_MS
    LSQ_groups{ii1}=find(groups==ii1);
end

end

