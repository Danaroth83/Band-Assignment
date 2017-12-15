function [Bands_to_sharpen,Bands_to_display]=load_Bands_to_sharpen(sensor,label,im_tag)

if nargin<=1, label='all'; end
if nargin<=2, im_tag=[]; end

if strcmpi(sensor,'Hyperion'), sensor='HYP'; end
if strcmpi(sensor,'CHRIS'), sensor='CHR'; end
if strcmpi(sensor,'IKO'), sensor='IKONOS'; end
if strcmpi(sensor,'GE1'), sensor='GeoEye1'; end
if strncmpi(im_tag,'Beijing',7) && strcmpi(sensor,'WV3'), sensor='WV34bands'; end

switch sensor
    case 'HYP'
        %%% Band to display options
        Bands_to_display=[16,23,29]; % Wavelength 508, 579, 640 (suggested on user guide)
        % Bands_to_display=[16,23,28]; % Wavelength 508, 579, 630 (minimum SAM test with MS)
        % Bands_to_display=[16,23,32]; % Wavelength 508, 579, 671 (best match to MS-ALI)
        % Bands_to_display=[13,18,30]; % Wavelength 478, 529, 651 (middle of color spectrum)       
        % Bands_to_display=[16,22,27]; % Wavelength 508, 569, 620 (original toolbox)
        % Bands_to_display=[14,18,33]; % Wavelength 488, 529, 681 (lowest visual fidelty within overlap bands)
        % Bands_to_display=[14,20,32]; % Wavelength 488, 548, 671
        if strcmpi(label,'all')
            Bands_to_sharpen=1:242;
        elseif strcmpi(label,'unnoise')
            if strncmpi(im_tag,'Beijing',7)
                Bands_to_sharpen=[9:57,82:97,99:119,134:164,182:184,187:220];
            elseif strncmpi(im_tag,'Sydney',6)
                Bands_to_sharpen=[9:57,82:97,99:119,134:164,182:184,187:220]; %Unchecked
            elseif strncmpi(im_tag,'RdJ',3)
                Bands_to_sharpen=[9:57,82:97,99:119,134:164,182:184,187:220]; %Unchecked
            elseif strncmpi(im_tag,'SaoPaulo',8)
                Bands_to_sharpen=[9:57,82:96,101:118,134:164,183,187:220];
            elseif strncmpi(im_tag,'SanFrancisco',12)
                Bands_to_sharpen=[11:57,82:96,101:119,134:164,183,188:189,193:215,217];
            elseif strncmpi(im_tag,'Sofia',5)
                Bands_to_sharpen=[11:57,82:96,101:119,134:164,183:184,187:220];
            elseif strncmpi(im_tag,'Sudbury',7)
                Bands_to_sharpen=[11:57,82:97,99:119,133:164,188:189,193:218];
            else
                Bands_to_sharpen=[9:57,82:97,99:119,134:164,182:184,187:220];
            end
        elseif strcmpi(label,'VNIR')
            Bands_to_sharpen=load_Bands_to_sharpen(sensor,'unnoise',im_tag);
            Bands_to_sharpen=intersect(Bands_to_sharpen,9:57);
        elseif strcmpi(label,'16')
            Bands_to_sharpen=[14:16,18:24,28:33];
        elseif strcmpi(label,'MSovlp')
            [~,sensor_ovl]=load_default_sensor(place);
            Bands_to_sharpen=load_Band_overlap(sensor,sensor_ovl,'HS','MS');
            Bands_to_sharpen=intersect(Bands_to_sharpen,load_Bands_to_sharpen(sensor,'16',im_tag));
        elseif strcmpi(label,'MSovlp2')
            Bands_to_sharpen=load_Band_overlap(sensor,'ALI','HS','MS');
            Bands_to_sharpen=intersect(Bands_to_sharpen,load_Bands_to_sharpen(sensor,'16',im_tag));
        elseif strcmpi(label,'MSovlpVNIR')
            [~,sensor_ovl]=load_default_sensor(place);
            Bands_to_sharpen=load_Band_overlap(sensor,sensor_ovl,'HS','MS');
            Bands_to_sharpen=intersect(Bands_to_sharpen,load_Bands_to_sharpen(sensor,'VNIR',im_tag));
        elseif strcmpi(label,'MSovlp2VNIR')
            Bands_to_sharpen=load_Band_overlap(sensor,'ALI','HS','MS');
            Bands_to_sharpen=intersect(Bands_to_sharpen,load_Bands_to_sharpen(sensor,'VNIR',im_tag));
        elseif strcmpi(label,'RGB')
            Bands_to_sharpen=[16,23,29];
        else
            Bands_to_sharpen=14:33;
        end
    case 'ALI'
        Bands_to_display=2:4;
        if any(strcmpi(label,{'all','unnoise'}))
            Bands_to_sharpen=1:9;
        elseif strcmpi(label,'VNIR')
            Bands_to_sharpen=1:6;
        elseif strcmp(label,'16')
            Bands_to_sharpen=1:4;
        elseif strcmp(label,'RGB')
            Bands_to_sharpen=2:4;
        else
            Bands_to_sharpen=1:6;
        end
    case 'WV2'
        if strcmpi(label,'RGB')
            Bands_to_sharpen=[1,3,5];
        else
            Bands_to_sharpen=1:8;
        end
        Bands_to_display=[1,3,5];
    case 'WV3'
        if strncmpi(im_tag,'Sydney',6) || strncmpi(im_tag,'RdJ',3)
            if strcmpi(label,'RGB')
                Bands_to_sharpen=[2,3,5];
            else
                Bands_to_sharpen=1:8;
            end
            Bands_to_display=[2,3,5];
        else
            if strcmpi(label,'RGB')
                Bands_to_sharpen=1:3;
            else
                Bands_to_sharpen=1:4;
            end
            Bands_to_display=1:3;
        end
    case 'ROSIS'
        if strcmpi(label,'32')
            Bands_to_sharpen=12:2:74; %32
            Bands_to_display=[16,30,60]; % Wavelength 491, 547, 667
        elseif strcmpi(label,'8')
            Bands_to_sharpen=10:10:80; %8
            Bands_to_display=[20,30,60];
        elseif strcmpi(label,'16')
            Bands_to_sharpen=10:5:85; %16
            Bands_to_display=[15,30,60];
        elseif strcmpi(label,'last')
            Bands_to_sharpen=51:1:66; %16
            Bands_to_display=[51,56,60];
        elseif strcmpi(label,'RGB')
            Bands_to_sharpen=[16,30,60];
            Bands_to_display=[16,30,60]; % Wavelength 491, 547, 667
        else
            Bands_to_sharpen=1:103;
            Bands_to_display=[16,30,60];
        end
    case 'AVIRIS'
        Bands_to_display=[3,9,21]; % Bands ?? watch out: to choose
        if strcmpi(label,'all')
            Bands_to_sharpen=1:176; % All bands
        elseif strcmpi(label,'4')
            Bands_to_sharpen=[3,9,15,21];
        elseif strcmpi(label,'ovlp')
            Bands_to_sharpen=1:41;
        elseif strcmpi(label,'RGB')
            Bands_to_sharpen=[3,9,21];
        else
            Bands_to_sharpen=1:2:32;
        end
    case 'CHR'
        if strcmpi(label,'16')
            Bands_to_sharpen=3:17; 
        elseif strcmpi(label,'RGB')
            Bands_to_sharpen=[7,9,11];
        else
            Bands_to_sharpen=1:18;
        end
        Bands_to_display=[7,9,11];
    case 'DE2'
        Bands_to_sharpen=1:4;
        Bands_to_display=[4,3,2];
    otherwise
        if strcmpi(label,'RGB')
           Bands_to_sharpen=1:3;
        else
            Bands_to_sharpen=1:4;
        end
        Bands_to_display=1:3;
end
        
end
            