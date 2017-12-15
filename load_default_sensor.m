function [sensor_lq,sensor_hq]=load_default_sensor(place)

if strncmpi(place,'RdJ',3) || strncmpi(place,'Sydney',6)
    sensor_lq='HYP';
    sensor_hq='WV3';
elseif strncmpi(place,'Beijing',7)
    sensor_lq='HYP';
    sensor_hq='WV34bands';
elseif strncmpi(place,'SanFrancisco',12)
    sensor_lq='HYP';
    sensor_hq='QB';
elseif strncmpi(place,'SaoPaulo',8)
    sensor_lq='HYP';
    sensor_hq='IKONOS';
elseif strncmpi(place,'Sofia',5) || strncmpi(place,'Sudbury',7)
    sensor_lq='HYP';
    sensor_hq='ALI';    
elseif strncmpi(place,'WV2',3) || strncmpi(place,'Rome',4) || strncmpi(place,'Rio',3)
    sensor_lq='WV2';
    sensor_hq='WV2';
elseif strncmpi(place,'China',5)
    sensor_lq='IKONOS';
    sensor_hq='IKONOS';
elseif strncmpi(place,'Hobart',6)
    sensor_lq='GeoEye1';
    sensor_hq='GeoEye1';
elseif strncmpi(place,'Indianapolis',12)
    sensor_lq='QB';
    sensor_hq='QB';
elseif strncmpi(place,'Tls',3) || strncmpi(place,'Toulouse',8)
    sensor_lq='IKONOS';
    sensor_hq='IKONOS';
elseif strncmpi(place,'CHR',3)
    sensor_lq='CHR';
    sensor_hq='CHR';
elseif strncmpi(place,'HYP',3) || strncmpi(place,'Paris',5)
    sensor_lq='HYP';
    sensor_hq='ALI';
elseif strncmpi(place,'ROSIS',5)
    sensor_lq='ROSIS';
    sensor_hq='ROSIS';
elseif strncmpi(place,'Moffett',7)
    sensor_lq='AVIRIS';
    sensor_hq='AVIRIS';
elseif strncmpi(place,'Vancouver',9)
    sensor_lq='DE2';
    sensor_hq='DE2';
end