clearvars; close all;
% im_tag='Beijing_cut4_HYP_ALI_ALI_VNIR';
im_tag='SanFrancisco_cut3_HYP_ALI_ALI_VNIR';
% im_tag='Sofia_cut2_HYP_ALI_ALI_VNIR';
flag_runalgorithms=1;

currentFolder=pwd;
cd ..; executiveFolder=pwd;

if flag_runalgorithms==1
    grouping_list={'SAMHS','LSQ','CC','MSD','CEN'};
    methods_fus={'EXP','Brovey','G-BDSD','GSA','ATWT','MTF-GLP-CBD','MTF-GLP-HPM','CNMF','HySure'};
    methods_opt={'ATWT','MTF-GLP-HPM','MTF-GLP-CBD'};
    qindex={'q2n','SAM','ERGAS'};
    
    MR=[]; MI=[];
    SCC_b_out=[]; ERGAS_b_out=[]; Q_b_out=[]; Band_assignment_out=[];
    time_grouping_mat=zeros(1,numel(grouping_list));
    time_mat=zeros(numel(methods_fus),numel(grouping_list));
    for ii=1:numel(grouping_list)
        [MI_temp,MR_temp,tit1,col1]=start_RR(im_tag,'M','ratio',3,'ratio_tgt_MS',3,'grouping',grouping_list{ii},...
            'interpolation',0,'methods',methods_fus,'val',qindex,'printMAT',1,'Qblocks_size',32,'draw',1);
        MI=cat(5,MI,MI_temp);
        MR=cat(3,MR,MR_temp);
        load('Output.mat');
        time_grouping_mat(ii)=time_grouping;
        time_mat(1:size(MI,4),ii)=time; time_mat(size(MI,4)+1:end,:)=[];
        SCC_b_out=cat(3,SCC_b_out,SCC_b);
        ERGAS_b_out=cat(3,ERGAS_b_out,ERGAS_b);
        Q_b_out=cat(3,Q_b_out,Q_b);
        Band_assignment_out=cat(1,Band_assignment_out,zeros(1,length(Bands_to_sharpen)));
        for i2=1:numel(Band_grouped), Band_assignment_out(ii,Band_grouped{i2})=Bands_to_sharpen_MS(i2); end
        close all;
    end
    %[opt_Q,opt_ERGAS,opt_SCC,opt_grouping_Q,opt_grouping_ERGAS,opt_grouping_SCC]=start_optMulti_RR('Beijing_ALI_cut2_VNIR','ratio',3,'ratio_tgt_MS',3,'methods',methods_opt);
    filename_output=[currentFolder,'\Output\',im_tag];
    save([filename_output,'.mat'],'grouping_list','methods_fus','methods_opt','MR','MI','tit1','col1','SCC_b_out','ERGAS_b_out','Q_b_out','Band_assignment_out',...
    'sensor_PAN','ratio','Band_overlap_MS','Bands_to_sharpen','I_GT','time_mat','time_grouping_mat');
else
    filename_output=[currentFolder,'\Output\',im_tag];
    load([filename_output,'.mat']);
end

idx=1;
tit=cell(size(MR,1)*size(MR,3),1);
MR_out=zeros(size(MR,1)*size(MR,3),size(MR,2));
cell_grouping=zeros(1,numel(tit1));
for ii=1:numel(tit1)
    if any(strcmpi(tit1{ii},{'EXP','CNMF','HySure'}))
        MR_out(idx,:)=mean(MR(ii,:,:),3);
        tit{idx}=sprintf('\\\\textbf{%s}&',tit1{ii});
        cell_grouping(ii)=1;
        idx=idx+1;
    else
        MR_out(idx:idx+numel(grouping_list)-1,:)=squeeze(MR(ii,:,:)).';
        tit{idx}=sprintf('\\\\multirow{%d}{*}{\\\\textbf{%s}}&%s',numel(grouping_list),tit1{ii},grouping_list{1});
        tit(idx+1:idx+numel(grouping_list)-1)=cellfun(@(x) sprintf('&%s',x),grouping_list(2:end),'Un',0);
        cell_grouping(ii)=numel(grouping_list);
        idx=idx+numel(grouping_list);
    end
end
MR_out=MR_out(1:idx-1,:);
tit=tit(1:idx-1);
col=col1;
col{1}=sprintf('&%s',col1{1});

fid=fopen([filename_output,'.tex'],'w');
fclose(fid);
matrix2latex_mod3(MR_out,[filename_output,'.tex'],'rowLabels',tit,...
    'columnLabels',col,'alignment','c','format','%.4f',...
    'highlight',matrix2latex_highlightoption(col),'gRow',cell_grouping,...
    'hRow',cell_grouping,'hCol',ones(1,length(col)));

% save([filename_output,'_v2.mat'],'tit','col','MR_out','SCC_b_out','ERGAS_b_out','Q_b_out','Band_assignment_out',...
%    'opt_Q','opt_ERGAS','opt_SCC','opt_grouping_Q','opt_grouping_ERGAS','opt_grouping_SCC',...
%     'sensor_PAN','ratio','Band_overlap_MS','Bands_to_sharpen','methods_fus','methods_opt');

save([filename_output,'_v2.mat'],'tit','col','MR_out','SCC_b_out','ERGAS_b_out','Q_b_out','Band_assignment_out',...
'sensor_PAN','ratio','Band_overlap_MS','Bands_to_sharpen','methods_fus','methods_opt');

load([filename_output,'.mat']);

BtS=[37,26,17];
method_vis='GSA';
idx_EXP=find(strcmpi(methods_fus,'EXP'));
idx_met=find(strcmpi(methods_fus,method_vis));

I_MS=MI(:,:,:,idx_EXP,1);
I_Fus=squeeze(MI(:,:,:,idx_met,:));
I_Det=cat(4,I_GT,I_Fus)-repmat(I_MS,[1,1,1,size(I_Fus,4)+1]);

BtS=arrayfun(@(x) find(Bands_to_sharpen==x),BtS,'Un',1);

grouping_list_2=[{'GT'},grouping_list];

cd(executiveFolder); cd Visualization;
I_Det_out=showImagesAll(I_Det,grouping_list_2,BtS,0,[],0);
for ii=1:size(I_Det_out,4)
    printLocal(I_Det_out(:,:,:,ii),[filename_output,'_',grouping_list_2{ii},'_',method_vis]);
end


cd(executiveFolder);
choice_group=1:numel(grouping_list);
Matrix_times=[time_grouping_mat(choice_group); time_mat(:,choice_group)];
Matrix_times=[Matrix_times,mean(Matrix_times,2)];
time_columnlabel=[grouping_list(choice_group),'Mean'];
time_rowlabel=['Grouping',tit1]';
fid=fopen([filename_output,'_time.tex'],'w');
fclose(fid);
matrix2latex_mod3(Matrix_times,[filename_output,'_time.tex'],'rowLabels',time_rowlabel,...
    'columnLabels',time_columnlabel,'alignment','c','format','%.4f');
save([filename_output,'_time.mat'],'Matrix_times','time_rowlabel','time_columnlabel');
cd(currentFolder);