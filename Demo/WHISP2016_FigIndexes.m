function [band_assignment,Q_choice,ERGAS_choice,SCC_choice]=...
        WHISP2016_FigIndexes(place,choices,colors,method)
    
%     function [band_assignment_Q_choice,...
%     %,band_assignment_ERGAS_choice,band_assignment_SCC_choice,
%     Q_choice,ERGAS_choice,SCC_choice]=...
%         WHISP2016_FigIndexes(place,choices,colors,method,sensor_inp,ratio_inp)
% Example
% [a,b,c,d]=WHISP2016_FigIndexes('Beijing',{'SCCFR','SAMHSFR','HySh','optQ','closest'},{'b-','r-','g-','k--','c--'},'MTF-GLP-CBD','Not_Ali',3);
% % [bg1,bg2]=Show_SAMclassifierdata(place,choice1,choice2,method,sensor,ratio)
% % 
% % ove choice1 e choice2 possono avere i seguenti valori:
% % 'SAMf': La versione originale dell'algoritmo SAM-AA
% % 'SAMs': La versione snella del SAM-AA
% % 'SAMsFR': Versione snella del SAM-AA fatta a FR
% % 'SAMHS': La versione proposta da voi per il SAM-AA
% % 'SAMHSFR': Come sopra ma a FR
% % 'CC': L'algoritmo CC-AA di Aiazzi
% % 'CCFR': Lo stesso a FR
% % 'closest': L'assegnazione closest band
% % 'HySh': L'algoritmo di hypersharpening che simula una PAN con regressione lineare proposto da Aiazzi
% % 'opt': E' l'ottimo per ognuno degli indici
% % 'optQ': Raggruppamento ottimo per l'UIQI
% % 'optERGAS': Raggruppamento ottimo per l'ERGAS
% % 'optSCC': Raggruppamento ottimo per l'SCC

% if nargin<=5, ratio_inp=3; end
% if nargin<=4, sensor_inp='notALI'; end
% if nargin<=3, method='MTF-GLP-HPM'; end
% if nargin<=2, choice2='SAM-AA'; end
% if nargin<=1, choice1='CC-AA'; end
% if nargin<=0, place='Beijing'; end
% keyboard
% close all;

method=strrep(method,'_','-');
load(['Output/',place,'_v2.mat']);
method_opt=0;

if ~ischar(method), error('Method should be a string'); end
flag_done=0;
for ii=1:numel(methods_opt)
    if strcmpi(methods_opt{ii},method), method_opt=ii; end
end
for ii=1:numel(methods_fus)
    if strcmpi(methods_fus{ii},method)
        method_fus=ii;
        flag_done=1;
    end
end
if flag_done==0, error('Unsupported method'); end


band_assignment=zeros(numel(choices),size(Band_assignment_out,2));
Q_choice=zeros(numel(choices),size(Band_assignment_out,2));
ERGAS_choice=zeros(numel(choices),size(Band_assignment_out,2));
SCC_choice=zeros(numel(choices),size(Band_assignment_out,2));
legend_choice=cell(1,numel(choices));
        
for ii=1:numel(choices)
    if any(strcmpi(choices{ii},{'Optimum','opt','optQ','optERGAS','optSCC'}))
        if method_opt==0, error('Optimum is only available for MRA methods'); end
        if strcmpi(choices{ii},'optERGAS')
            band_assignment(ii,:)=opt_grouping_ERGAS(method_opt,:);
        elseif strcmpi(choices{ii},'optSCC')
            band_assignment(ii,:)=opt_grouping_SCC(method_opt,:);
        else
            band_assignment(ii,:)=opt_grouping_Q(method_opt,:);
        end
        Q_choice(ii,:)=opt_Q(method_opt,:);
        ERGAS_choice(ii,:)=opt_ERGAS(method_opt,:);
        SCC_choice(ii,:)=opt_SCC(method_opt,:);
        legend_choice{ii}='Optimum';
    else
        if any(strcmpi(choices{ii},{'CC-AA','CC','SCC-AA','SCC','CCRR-AA','CCRR'}))
            choices{ii}='CC';
            legend_choice{ii}='CC-AA (RR)';
        elseif any(strcmpi(choices{ii},{'SAMHS-AA','SAMHS','SAMHSRR-AA','SAMHSRR','SAM-AA','SAM'}))
            choices{ii}='SAMHS';
            legend_choice{ii}='SAMHS-AA (RR)';
        elseif any(strcmpi(choices{ii},{'Closest','BAD','MSD','Bands'}))
            choices{ii}='MSD';
            legend_choice{ii}='MSD';
        elseif any(strcmpi(choices{ii},{'LSQ','MMSE'}))
            choices{ii}='LSQ';
            legend_choice{ii}='MMSE';
        elseif any(strcmpi(choices{ii},{'CEN','Centroid'}))
            choices{ii}='CEN';
            legend_choice{ii}='Centroid';
        end
        group=0;
        group=find(strcmpi(choices{ii},grouping_list));
        if group==0, error('Assignment method is not available; please check name'); end
        band_assignment(ii,:)=Band_assignment_out(group,:);
        Q_choice(ii,:)=Q_b_out(:,method_fus,group).';
        ERGAS_choice(ii,:)=ERGAS_b_out(:,method_fus,group).';
        SCC_choice(ii,:)=SCC_b_out(:,method_fus,group).';
    end

end


% keyboard

color_bands={[0,0,1],[0,1,0],[1,0,0],[1,0,1],[0,1,1],[1,1,0],[1,0.5,0],[0,0.5,1]};
Transparency=0.1;
% Region_overlap=unique(cell2mat(band_overlap_cell));
% color_bands=['b','g','r','m', 'c','y'];
% RGB=[0 0 1;0 1 0;1 0 0;1 0 1;0 1 1;1 1 0;1 0.5 0;0 0.5 1];
%  HSV=rgb2hsv(RGB);
% HSV(:,2)=0.2*HSV(:,2);
% RGB2=hsv2rgb(HSV);
%  color_bands= mat2cell(RGB2,ones(1,size(RGB2,1)), 3)';
%  keyboard
%      colour_teal = [18 150 155] ./ 255;
% colour_lightgreen = [94 250 81] ./ 255;
% colour_green = [12 195 82] ./ 255;
% colour_lightblue = [8 180 238] ./ 255;
% colour_darkblue = [1 17 181] ./ 255;
% colour_yellow = [251 250 48] ./ 255;
% colour_peach = [251 111 66] ./ 255;
%  color_bands={colour_lightblue,colour_lightgreen,[1 0 0],[1 0 1],[0 1 1],[1 1 0],[1 0.5 0],[0 0.5 1]}; 
band_overlap_cell=cellfun(@(x) Bands_to_sharpen(x), Band_overlap_MS,'Un',0);

hq=figure; 
 hold on;
for ii=1:1:numel(choices),
    plot(Bands_to_sharpen,Q_choice(ii,:),colors{ii},'LineWidth',1.5);
end
legend(legend_choice,'Location','Southwest')
title(sprintf(['Sensor: ',sensor_PAN,'; Method: ',method,'; Ratio: %d'],ratio));
xlabel('Bands');
ylabel('Q');
Ylimits = get(gca,'YLim');
for jj=1:length(band_overlap_cell)
    Region_overlap=cell2mat(band_overlap_cell(jj));
    for ii=1:length(Region_overlap)
        patch('Faces',[1 2 3 4],'Vertices',[Region_overlap(ii)-0.5 Ylimits(1);Region_overlap(ii)+0.5 Ylimits(1);Region_overlap(ii)+0.5 Ylimits(2);Region_overlap(ii)-0.5 Ylimits(2)],'FaceColor',color_bands{jj},'EdgeColor','none','FaceAlpha',Transparency)
    end
end
ylim(Ylimits);
% for ii1=1:length(band_assignment_Q_choice1)
%     text(Bands_to_sharpen(ii1),Ylimits(2)-0.025*(Ylimits(2)-Ylimits(1)),num2str(band_assignment_Q_choice1(ii1)),'Color','red','FontSize',6,'HorizontalAlignment','center');
%     text(Bands_to_sharpen(ii1),Ylimits(1)+0.025*(Ylimits(2)-Ylimits(1)),num2str(band_assignment_Q_choice2(ii1)),'Color','blue','FontSize',6,'HorizontalAlignment','center');
% end


he=figure; 
 hold on;
for ii=1:1:numel(choices),
    plot(Bands_to_sharpen,ERGAS_choice(ii,:),colors{ii},'LineWidth',1.5);
end
legend(legend_choice,'Location','Southwest')
title(sprintf(['Sensor: ',sensor_PAN,'; Method: ',method,'; Ratio: %d'],ratio));
xlabel('Bands');
ylabel('ERGAS');
Ylimits = get(gca,'YLim');
% for ii=1:length(Region_overlap)
%     patch('Faces',[1 2 3 4],'Vertices',[Region_overlap(ii)-0.5 Ylimits(1);Region_overlap(ii)+0.5 Ylimits(1);Region_overlap(ii)+0.5 Ylimits(2);Region_overlap(ii)-0.5 Ylimits(2)],'FaceColor','g','EdgeColor','none','FaceAlpha',0.3)
% end
for jj=1:length(band_overlap_cell)
    Region_overlap=cell2mat(band_overlap_cell(jj));
    for ii=1:length(Region_overlap)
        patch('Faces',[1 2 3 4],'Vertices',[Region_overlap(ii)-0.5 Ylimits(1);Region_overlap(ii)+0.5 Ylimits(1);Region_overlap(ii)+0.5 Ylimits(2);Region_overlap(ii)-0.5 Ylimits(2)],'FaceColor',color_bands{jj},'EdgeColor','none','FaceAlpha',Transparency)
    end
end
ylim(Ylimits);
% for ii1=1:length(band_assignment_Q_choice1)
%     text(Bands_to_sharpen(ii1),Ylimits(2)-0.025*(Ylimits(2)-Ylimits(1)),num2str(band_assignment_Q_choice1(ii1)),'Color','red','FontSize',6,'HorizontalAlignment','center');
%     text(Bands_to_sharpen(ii1),Ylimits(1)+0.025*(Ylimits(2)-Ylimits(1)),num2str(band_assignment_Q_choice2(ii1)),'Color','blue','FontSize',6,'HorizontalAlignment','center');
% end


% hs=figure; 
%  hold on;
% for ii=1:1:numel(choices),
%     plot(Bands_to_sharpen,SCC_choice(ii,:),colors{ii},'LineWidth',1.5);
% end
% legend(legend_choice,'Location','Southwest')
% title(sprintf(['Sensor: ',sensor,'; Method: ',titleImages{method},'; Ratio: %d'],ratio_inp));
% xlabel('Bands');
% ylabel('SCC');
% Ylimits = get(gca,'YLim');
% % for ii=1:length(Region_overlap)
% %     patch('Faces',[1 2 3 4],'Vertices',[Region_overlap(ii)-0.5 Ylimits(1);Region_overlap(ii)+0.5 Ylimits(1);Region_overlap(ii)+0.5 Ylimits(2);Region_overlap(ii)-0.5 Ylimits(2)],'FaceColor','g','EdgeColor','none','FaceAlpha',0.3)
% % end
% for jj=1:length(band_overlap_cell)
%     Region_overlap=cell2mat(band_overlap_cell(jj));
%     for ii=1:length(Region_overlap)
%         patch('Faces',[1 2 3 4],'Vertices',[Region_overlap(ii)-0.5 Ylimits(1);Region_overlap(ii)+0.5 Ylimits(1);Region_overlap(ii)+0.5 Ylimits(2);Region_overlap(ii)-0.5 Ylimits(2)],'FaceColor',color_bands{jj},'EdgeColor','none','FaceAlpha',Transparency)
%     end
% end
% ylim(Ylimits);
% % for ii1=1:length(band_assignment_Q_choice1)
% %     text(Bands_to_sharpen(ii1),Ylimits(2)-0.025*(Ylimits(2)-Ylimits(1)),num2str(band_assignment_Q_choice1(ii1)),'Color','red','FontSize',6,'HorizontalAlignment','center');
% %     text(Bands_to_sharpen(ii1),Ylimits(1)+0.025*(Ylimits(2)-Ylimits(1)),num2str(band_assignment_Q_choice2(ii1)),'Color','blue','FontSize',6,'HorizontalAlignment','center');
% % end
% keyboard

hb=figure; 
 hold on;
for ii=1:1:numel(choices),
    plot(Bands_to_sharpen, band_assignment(ii,:),colors{ii},'LineWidth',1.5);
end
legend(legend_choice,'Location','Northwest')
title(sprintf(['Sensor: ',sensor_PAN,'; Method: ',method,'; Ratio: %d'],ratio));
xlabel('Bands');
ylabel('Assignment');
% Ylimits = get(gca,'YLim');
Ylimits =([min(band_assignment(:))-1 max(band_assignment(:))+1]);

% for ii=1:length(Region_overlap)
%     patch('Faces',[1 2 3 4],'Vertices',[Region_overlap(ii)-0.5 Ylimits(1);Region_overlap(ii)+0.5 Ylimits(1);Region_overlap(ii)+0.5 Ylimits(2);Region_overlap(ii)-0.5 Ylimits(2)],'FaceColor','g','EdgeColor','none','FaceAlpha',0.3)
% end
for jj=1:length(band_overlap_cell)
    Region_overlap=cell2mat(band_overlap_cell(jj));
    for ii=1:length(Region_overlap)
        patch('Faces',[1 2 3 4],'Vertices',[Region_overlap(ii)-0.5 Ylimits(1);Region_overlap(ii)+0.5 Ylimits(1);Region_overlap(ii)+0.5 Ylimits(2);Region_overlap(ii)-0.5 Ylimits(2)],'FaceColor',color_bands{jj},'EdgeColor','none','FaceAlpha',Transparency)
    end
end
ylim(Ylimits);
% for ii1=1:length(band_assignment_Q_choice1)
%     text(Bands_to_sharpen(ii1),Ylimits(2)-0.025*(Ylimits(2)-Ylimits(1)),num2str(band_assignment_Q_choice1(ii1)),'Color','red','FontSize',6,'HorizontalAlignment','center');
%     text(Bands_to_sharpen(ii1),Ylimits(1)+0.025*(Ylimits(2)-Ylimits(1)),num2str(band_assignment_Q_choice2(ii1)),'Color','blue','FontSize',6,'HorizontalAlignment','center');
% end
