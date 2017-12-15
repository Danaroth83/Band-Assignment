function [ spatial_res,radiometric_res ] = load_resolution( sensor, im_tag, type )
%LOAD_RESOLUTION Loads the spatial resolution in meters of specific satellite sensors
if nargin<=2 || isempty(type), type='MS'; end
if nargin<=1, im_tag=[]; end
if strcmpi(sensor,'IKO'), sensor='IKONOS'; end
if strcmpi(sensor,'Hyperion'), sensor='HYP'; end
if strcmpi(sensor,'CHRIS'), sensor='CHR'; end
if strcmpi(sensor,'GE1'), sensor='GeoEye1'; end

if isequal(type,1) || strcmp(type,'PAN')
    switch sensor
        case 'QB'
            spatial_res=0.6;
        case 'IKONOS'
            spatial_res=0.8;
        case 'WV1'
            spatial_res=0.5;
        case 'WV2'
            spatial_res=0.4;
        case {'WV3','WV34bands'}
            if strncmpi(im_tag,'Adelaide',8)
                spatial_res=0.3;
            elseif strncmpi(im_tag,'Beijing',7)
                spatial_res=0.4;
            elseif strncmpi(im_tag,'Sydney',6)
                spatial_res=0.5;
            elseif strncmpi(im_tag,'RdJ',3)
                spatial_res=0.3;
            else
                spatial_res=0.4;
            end
        case 'GeoEye1'
            spatial_res=0.5;
        case 'DE2'
            spatial_res=1;
        case {'ALI','HYP'}
            spatial_res=10;
        otherwise
            spatial_res=NaN;
    end
else
    switch sensor
        case 'QB'
            spatial_res=2.4;
        case 'IKONOS'
            spatial_res=3.2;
        case 'WV1'
            spatial_res=2;
        case 'WV2'
            spatial_res=1.6;
        case {'WV3','WV34bands'}
            if strncmpi(im_tag,'Adelaide',8)
                spatial_res=1.2;
            elseif strncmpi(im_tag,'Beijing',7)
                spatial_res=1.6;
            elseif strncmpi(im_tag,'Sydney',6)
                spatial_res=2;
            elseif strncmpi(im_tag,'RdJ',3)
                spatial_res=1.2;
            else
                spatial_res=1.6;
            end
        case 'GeoEye1'
            spatial_res=2;
        case 'DE2'
            spatial_res=4;
        case 'ALI'
            spatial_res=30;
        case 'HYP'
            spatial_res=30;
        case 'AVIRIS'
            spatial_res=30;
        otherwise
            spatial_res=NaN;
    end
end

if nargout>=2
    switch sensor
        case {'WV1','WV2','WV3','WV34bands','IKONOS','QB','GeoEye1','DE2'}
            radiometric_res=11;
        case {'ALI','HYP','ROSIS'}
            radiometric_res=15;
        case 'CHR'
            radiometric_res=16;
        case 'AVIRIS'
            radiometric_res=13;
        otherwise
            radiometric_res=NaN;
    end
end
end

