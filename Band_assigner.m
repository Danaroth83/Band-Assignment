function [Groups,time_grouping]=Band_assigner(flag_grouping,I_PAN,I_MS_LR,ratio,I_MS,GNyq_MS,band_overlap_cell,Bands_to_sharpen,Bands_to_sharpen_MS,sensor_HS,sensor_MS,im_tag)

cd Assignment
t2=tic;
if isequal(flag_grouping,2) || strcmpi(flag_grouping,'CC')
    Groups=CC_classifier_v2(I_PAN,I_MS_LR,ratio,GNyq_MS);
elseif isequal(flag_grouping,3) || strcmpi(flag_grouping,'MSD')
    Groups=BaDi_classifier(sensor_HS,sensor_MS,Bands_to_sharpen,Bands_to_sharpen_MS,im_tag);
elseif isequal(flag_grouping,5) || strcmpi(flag_grouping,'SAMHS')
    Groups=SAMHS_classifier_v2(I_PAN,I_MS_LR,ratio,GNyq_MS);
elseif isequal(flag_grouping,20) || strcmpi(flag_grouping,'Centroid')
    Groups=Centroid_classifier(sensor_HS,sensor_MS,Bands_to_sharpen,Bands_to_sharpen_MS,im_tag);
elseif isequal(flag_grouping,26) || strcmpi(flag_grouping,'LSQ')
    Groups=LSQ_classifier(I_PAN,I_MS);
else
    Groups=band_overlap_cell;
end
time_grouping=toc(t2);
cd ..

end