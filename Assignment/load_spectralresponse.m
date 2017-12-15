function [ central,dispersion,wavelength_nm,response ] = load_spectralresponse( sensor,type,im_tag,Bands )
%LOAD_SPECTRALRESPONSE
%   Loads the centroid and the dispersion of an assigned spectral response
if nargin<=1, type='MS'; end
if strcmpi(sensor,'Hyperion'), sensor='HYP'; end
if strcmpi(sensor,'HYP') && strcmpi(type,'PAN'), sensor='ALI'; end
if nargin<=2, im_tag=[]; end
if strncmpi(im_tag,'Beijing',7) && strcmpi(sensor,'WV3'), sensor='WV34bands'; end
if strcmpi(sensor,'DE2'), [central,dispersion]=load_wavelength(Bands,sensor,im_tag,type); return; end

currentFolder=pwd;
cd('../Relative Spectral Responses/');
load([sensor,'_Spectral_Responses.mat']);
cd(currentFolder);
if ~strcmp(type,'PAN')
    if nargin<=3
        if ~any(strcmpi(sensor,{'HYP','AVIRIS'}))
            Bands=1:size(Spectral_Responses_Matrix,1)-1;
        else
            Bands=1:size(Spectral_Responses_Matrix,1);
        end
    end
    Nbands=length(Bands);
    Nwavelengths=length(wavelength_nm);
    response=Spectral_Responses_Matrix(Bands,:);
    summation=sum(response,2);
    wavelength_mat=repmat(wavelength_nm.',[Nbands,1]);
    central=(sum(wavelength_mat.*response,2)./summation);
    dispersion=sqrt(sum(response.*(wavelength_mat-repmat(central,[1,Nwavelengths])).^2,2)./summation);
    dispersion=dispersion';
    central=central';
else
    response=Spectral_Responses_Matrix(end,:);
    summation=sum(response,2);
    central=sum(wavelength_nm.'.*response)/summation;
    dispersion=sqrt(sum(response.*(wavelength_nm.'-repmat(central,[1,length(wavelength_nm)])).^2)/summation);
end
end

