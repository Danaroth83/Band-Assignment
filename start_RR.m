function [MatrixImage,MatrixResults,titleImages,columnLabels]=start_RR(im_tag,shortname,varargin)

% Start function for sharpening of low spatial quality images and Reduced
% Resolution validation
% Usage: [MI,MR]=start_RR(im_tag,type,''ratio'',3)
% im_tag is a label identifying the dataset (eg:
% Beijing_cut2_HYP_WV3_WV3_VNIR is Beijing second cut dataset with an
% Hyperion HS to sharpen with a WV3 MS and PAN, considering only VNIR
% bands)
% type is the sharpening method (eg: 'P' for classical pansharpening, 'M'
% for HS/MS sharpening
% Dynamic inputs are allowed:
% 'ratio' is the scale ratio between the HQ image and the image to sharpen
% 'ratio_REF' is the target ratio of the sharpened image
% 'methods' is the list of fusion methods (eg: {'EXP','GSA','MTF-GLP-CBD'})
% 'grouping' is the HS/MS assignment method (eg: 'SAMHS','CC','OVLP')
% 'interp' is the interpolation method (eg: 0=bicubic, 1=fir1 kernel, 2=Lagrange kernel)
% 'draw' if 0, it bypasses drawing figures
% 'printMAT' if 1, it saves processed data in Output.mat
% 'printEPS' if 1, it saves printed images in /Outputs/
% Output:
% MI: 4D Matrix of output results (length x width x bands x fusion methods)
% MR: 2D Matrix of validation results (indexes x fusion methods)

if nargin<1, error('Image tag was not inputted'); end
if nargin<=1, shortname='P'; end
if numel(varargin)==1, varargin=varargin{1}; end
if rem(numel(varargin),2)~=0
    fprintf(['Usage: start_RR(im_tag,''P'',''ratio'',ratio,''ratio_GT'',ratio_GT,',...
        '''methods'',methods);\n''P'' is for MS/PAN, ''M'' for HS/MS sharpening\n']);
    error('Wrong input format');
end

setup_flags;

for ii=1:2:numel(varargin)
    pname=varargin{ii};
    pval=varargin{ii+1};
    if strcmpi(pname,'ratio')                                           % Ratio between reference and MS-LR
        ratio=pval;
    elseif any(strcmpi(pname,{'grouping','assignment','flag_grouping'}))% Assignment method
        flag_grouping=pval;
        flag_grouping_list=[0,1,2,3,4,5,6,7,8,9,10,12,14,15,20,21,22,23,24,25,26,27];
        grouping_label_list={'OVLP','SAM','CC','MSD','SAM2','SAMHS','SAMHSEXP','HYSH','CCEXP','SAMHSEXP2','CCHM','CCFR','SAMFR','SAMHSFR','CEN','DIS','CEN2','DIS2','CCOV','RSR','LSQ','LSQ2'};
        if any(strcmpi(flag_grouping,grouping_label_list))
            flag_grouping=flag_grouping_list(strcmpi(flag_grouping,grouping_label_list));
        end
    elseif any(strcmpi(pname,{'ratio_tgt_PAN','ratio_PAN'}))            % Initial target scale ratio of PAN
        ratio_tgt_PAN=pval;
    elseif any(strcmpi(pname,{'ratio_tgt_MS','ratio_MS'}))              % Initial target scale ratio of MS
        ratio_tgt_MS=pval; 
    elseif any(strcmpi(pname,{'ratio_tgt_REF','ratio_REF','ratio_GT'})) % Ratio between reference and MS-LR
        ratio_tgt_REF=pval;
    elseif any(strcmpi(pname,{'methods','methods_list'}))               % List of fusion methods (in a cell)
        methods_list=load_methodslist(pval);
    elseif any(strcmpi(pname,{'bandselect','flag_bandselect_val','bandselect_val'})) % 2=Validation only on non-overlapping bands
        flag_bandselect_val=pval;
    elseif any(strcmpi(pname,{'cut_bounds','cutbounds'}))               % If 1, it cuts bounds in validation
        flag_cutbounds_val=pval;
    elseif any(strcmpi(pname,{'thvalues'}))                             % If 1, it thresholds output values
        flag_thvalues_val=pval;
    elseif any(strcmpi(pname,{'val','RR_val','flag_RR_val'}))           % Chooses the method for RR validation
        if iscell(pval) || pval==6
            flag_RR_val=6;
            qindex_list=pval;
        else
            flag_RR_val=pval;
        end
    elseif any(strcmpi(pname,{'flag_interpolation','interpolation','upsample','interp'})) % Selects the interpolation type
        flag_interpolation=pval;
    elseif any(strcmpi(pname,{'draw','flag_draw'}))                     % If 0, it bypassses image visualization
        flag_draw=pval;
    elseif any(strcmpi(pname,{'printMAT','flag_printMAT'}))             % If 1, it saves data on a .mat file
        printMAT=pval;
    elseif any(strcmpi(pname,{'printEPS','flag_printEPS'}))             % If 1, it saves figures on .eps files
        printEPS=pval;
    elseif any(strcmpi(pname,{'Qblocks_size'}))
        Qblocks_size_setup=pval;
        flag_noautoblock=1;
    end
end
clear('pval','pname','varargin','flag_grouping_list','grouping_label_list');

switch shortname
    case 'P'
        disp('Pansharpening (Reduced Resolution)');
        Load_Dataset_Pansharpening_RR;
        Fusion_Algorithms_Pansharpening;
        Validation_RR;
    case 'M'
        disp('HS/MS Fusion (Reduced Resolution)');
        scale_ref=1;
        PAN_align=2;
        Load_Dataset_Multimodal_RR;
        Fusion_Algorithms_Multisharpening;
        Validation_RR;
end
