function [ min_sr,max_sr ] = load_minmaxspectralresponse( sensor,type,im_tag,Bands )
%LOAD_MINMAXSPECTRALRESPONSE
%   Loads the edges of a spectral response
if nargin<=1, type='MS'; end
if strcmpi(sensor,'Hyperion'), sensor='HYP'; end
if strcmpi(sensor,'IKO'), sensor='IKONOS'; end
if strcmpi(sensor,'HYP') && strcmpi(type,'PAN'), sensor='ALI'; end
if nargin<=2, im_tag=[]; end
if strncmpi(im_tag,'Beijing',7) && strcmpi(sensor,'WV3'), sensor='WV34bands'; end
if strcmpi(sensor,'AVIRIS') && strcmpi(type,'PAN')
    min_sr=400; max_sr=773; disp('Simulated AVIRIS PAN sensor was employed');
    return;
end
if strcmpi(sensor,'DE2'),
    if nargin<=3, Bands=1:4; end
    [central,dispersion]=load_wavelength(Bands,sensor,im_tag,type);
    min_sr=central-dispersion/2;
    max_sr=central+dispersion/2;
    return;
end


stopband_amp=0.5;

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
    response=Spectral_Responses_Matrix(Bands,:);
else
    Bands=1;
    response=Spectral_Responses_Matrix(end,:);
end

Nbands=length(Bands);
Nwavelengths=length(wavelength_nm);
maxamplitude=max(response,[],2);
higher=response>stopband_amp*repmat(maxamplitude,[1,Nwavelengths]);
min_sr=zeros(1,Nbands);
max_sr=zeros(1,Nbands);
for ii=1:Nbands
    min_sr(ii)=wavelength_nm(find(higher(ii,:),1));
    max_sr(ii)=wavelength_nm(find(higher(ii,:),1,'last'));
end
    
end