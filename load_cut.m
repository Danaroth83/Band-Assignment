function [hcut,vcut,edge_cut_lq,edge_cut_hq,Qblocks_size]=load_cut(im_tag,cut_label,sensor_lq,sensor_hq,type_hq)
if nargin<=3, sensor_hq=sensor_lq; end
if nargin<=4, type_hq='PAN'; end
if strcmpi(sensor_lq,'IKO'), sensor_lq='IKONOS'; end
if strcmpi(sensor_hq,'IKO'), sensor_hq='IKONOS'; end
if strcmpi(sensor_lq,'Hyperion'), sensor_lq='HYP'; end
if strcmpi(sensor_hq,'Hyperion'), sensor_hq='HYP'; end
if strcmpi(sensor_lq,'CHRIS'), sensor_lq='CHR'; end
if strcmpi(sensor_hq,'CHRIS'), sensor_hq='CHR'; end
if strcmpi(sensor_lq,'GE1'), sensor_lq='GeoEye1'; end
if strcmpi(sensor_hq,'GE1'), sensor_hq='GeoEye1'; end
if strcmpi(sensor_lq,'Deimos-2'), sensor_lq='DE2'; end
if strcmpi(sensor_hq,'Deimos-2'), sensor_hq='DE2'; end

ratio=load_resolution(sensor_lq,im_tag,'MS')/load_resolution(sensor_hq,im_tag,type_hq);
edge_cut_lq=0;
edge_cut_hq=0;
Qblocks_size=32;

if strncmpi(im_tag,'Beijing',7)
    if any(strcmpi(sensor_lq,{'HYP','ALI'}))
        Qblocks_size=27;
        edge_cut_lq=63;
        switch cut_label
            case 'cut2'
                hcut=68+(1:96); vcut=20+(1:144);      %Full image
            case 'cut3'
                hcut=128+(1:36); vcut=8+(1:156);
            case 'cut4'
                hcut=1:126; vcut=1:288;
            case 'cut5'
                hcut=4+(1:156); vcut=4+(1:156);
            otherwise
                hcut=116+(1:48); vcut=68+(1:96);
        end
    elseif any(strcmpi(sensor_lq,{'WV3','WV34bands'}))
        Qblocks_size=32;
        edge_cut_lq=44;
        edge_cut_hq=44*ratio;
        hcut=1458+(1:512); vcut=1716+(1:512);  %Bird's nest cut
        hcut=hcut-edge_cut_lq; vcut=vcut-edge_cut_lq;
    end
    if strcmpi(sensor_hq,'ALI')
        edge_cut_hq=63*ratio;
    end
elseif strncmpi(im_tag,'RdJ',3)
    if any(strcmpi(sensor_lq,{'HYP','ALI'}))
        Qblocks_size=27;
        edge_cut_lq=69;
        switch cut_label
            case 'cut2'
                hcut=2+(1:132); vcut=1+(1:108);  % Full Image
            otherwise
                hcut=12+(1:120); vcut=1+(1:108);
        end
    elseif strcmpi(sensor_lq,'WV3')
        Qblocks_size=32;
        edge_cut_lq=44;
        edge_cut_hq=44*ratio;
        hcut=1424+(1:512); vcut=44+(1:512);
        hcut=hcut-edge_cut_lq; vcut=vcut-edge_cut_lq;
    end
    if strcmpi(sensor_hq,'ALI')
        edge_cut_hq=69*ratio;
    end
elseif strncmpi(im_tag,'SaoPaulo',8)
    if any(strcmpi(sensor_lq,{'HYP','ALI'}))
        Qblocks_size=27;
        edge_cut_lq=138;
        switch cut_label
            case 'cut2'
                hcut=40+(1:192); vcut=40+(1:192);
            otherwise
                hcut=88+(1:144); vcut=40+(1:192);
        end
    elseif strcmpi(sensor_lq,'IKONOS')
        Qblocks_size=32;
        edge_cut_lq=44;
        edge_cut_hq=44*ratio;
        hcut=800+(1:256); vcut=800+(1:256);
        hcut=hcut-edge_cut_lq; vcut=vcut-edge_cut_lq;
    end
    if strcmpi(sensor_hq,'ALI')
        edge_cut_hq=138*ratio;
    end
elseif strncmpi(im_tag,'SanFrancisco',8)
    if any(strcmpi(sensor_lq,{'HYP','ALI'}))
        Qblocks_size=27;
        edge_cut_lq=0;
        switch cut_label
            case 'cut2'
                hcut=29:160; vcut=29:160;
                % hcut=25:160; vcut=25:160;
            case 'cut3'
                hcut=1:132; vcut=1:294;
            otherwise
                hcut=40+(1:120); vcut=88+(1:72);
        end
    elseif strcmpi(sensor_lq,'QB')
        Qblocks_size=32;
        edge_cut_lq=44;
        edge_cut_hq=44*ratio;
        switch cut_label
            case 'cut2'
                hcut=1376+(1:600); vcut=872+(1:600);
            otherwise
                hcut=1504+(1:256); vcut=1000+(1:256);
        end
        hcut=hcut-edge_cut_lq; vcut=vcut-edge_cut_lq;
    end
    if strcmpi(sensor_hq,'ALI')
        edge_cut_hq=0;
    end
elseif strncmpi(im_tag,'Sofia',5) || strncmpi(im_tag,'Sudbury',7)
    Qblocks_size=30;
    edge_cut_lq=60;
    edge_cut_hq=60*ratio;
    switch cut_label
        case 'cut2'
            if strncmpi(im_tag,'Sofia1GST',9) || strncmpi(im_tag,'Sudbury1GST',11)
                hcut=22:177; vcut=97:288;
            else
                hcut=1:156; vcut=1:192;
            end
        otherwise
            % hcut=1:192; vcut=97:288;   % Square; RR power of 2; black corner
            hcut=49:144; vcut=1:384;   % Rectangular; RR power of 2; no black corner
            % hcut=19:174; vcut=97:288;  % Rectangular; RR not power of 2; no black corner
            % hcut=17:176; vcut=113:272; % Square; RR not power of 2; no black corner
    end
elseif strncmpi(im_tag,'China',5)
    switch cut_label
        case 'cut2'
            hcut=22:53; vcut=22:53;
        otherwise
            hcut=1:75; vcut=1:75;
    end
elseif strncmpi(im_tag,'Hobart',6)
    if strncmpi(im_tag,'Hobart1_2',9) || strncmpi(im_tag,'Hobart2',7) || strcmpi(cut_label,'cut2')
        hcut=1:512; vcut=1:512;
    else
        hcut=1:128; vcut=1:128;
    end
elseif strncmpi(im_tag,'Rio',3)
    hcut=1:64; vcut=1:64;
elseif strncmpi(im_tag,'Rome',4)
    hcut=1:300; vcut=1:300;
elseif strncmpi(im_tag,'Indianapolis',12)
    if strncmpi(im_tag,'Indianapolis2',13)
        hcut=1:128; vcut=1:128;
    else
        hcut=1:25; vcut=1:25;
    end
elseif strncmpi(im_tag,'Tls',3) || strncmpi(im_tag,'Toulouse',8)
    hcut=1:128; vcut=1:128;
elseif strncmpi(im_tag,'CHR',3)
    hcut=1:100; vcut=1:100;
elseif strncmpi(im_tag,'HYP',3) || strncmpi(im_tag,'Paris',5)
    switch cut_label
        case 'cut2'
            hcut=181:280; vcut=31:130;
        otherwise
            hcut=141:290; vcut=11:160;
    end
elseif strncmpi(im_tag,'ROSIS',5)
    hcut=20+(1:256); vcut=50+(1:512);
elseif strncmpi(im_tag,'Moffett',7)
    switch cut_label
        case 'cut2'
            hcut=1:36; vcut=1:78; % Almost full image
        otherwise
            hcut=2+(1:32); vcut=7+(1:64);
    end
elseif strncmpi(im_tag,'Vancouver',9)
    hcut=299+(1:500); vcut=511+(1:500);
end
    
end