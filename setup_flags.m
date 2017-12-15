%% Cut Final Image
flag_cutbounds_vis = 0;  % Removes edges in visualization
flag_cutbounds_val = 0;  % Removes edges in quality assessment
dim_cut = 11;            % Amount of edge removal
flag_thvalues = 0;       % If 1, it thresholds output images to resolution range
flag_thvalues_vis = 0;   % If 1, it thresholds images for visualization
flag_thvalues_val = 0;   % If 1, it thresholds images for validation

%% Print Options
printEPS=0;              % Prints EPS files for figures
printTEX=1;              % Prints results in a TEX file
printTIFF=0;             % Prints TIFF files for figures
printMAT=1;              % Prints workspace variables in a .mat file
flag_draw=1;             % if 0, it bypasses generating and saving figures

%% Equalization options
flag_degradation=0;      % 0= non-degraded PAN; 1=degraded with imresize+imresize; 2=with MTF-PAN; 3=with ATWT; 4=with Hamming windowed sinc; 5=with MTF-MS
flag_equalization=0;     % 0=classic; 1=with GT; 2=regression based; 3=PCA regression based

%% Resizing options
%Upsampling
flag_interpolation=2;    % Upsampling option: 1=bicubic (imresize), 0=with a kernel generated with fir1, 2=with a Lagrange polynomial kernel

%Downsampling
flag_imresize_PAN=1;     % Downsampling of PAN when necessary 1=with imresize, 0=with MTF-filtered downsampling, 2=with HS's MTF emulated filter resize
flag_imresize_MS=0;      % Downsampling of MS when necessary 1=with imresize, 0=with MTF-filtered downsampling, 2=with HS's MTF emulated filter resize
flag_imresize_PAN_sim=1; % Downsampling of PAN to simulate a different resolution input PAN image
flag_imresize_MS_sim=0;  % Downsampling of PAN to simulate a different resolution input MS image
flag_imresize_PAN_RR=1;  % Resizing PAN for RR validation 1=with imresize, 0=with MTF-filtered downsampling
flag_imresize_REF_RR=1;  % Resizing REF for RR validation 1=with imresize, 0=with MTF-filtered downsampling
flag_imresize_MSREF_RR=0;% Resizing MSREF for RR validation 1=with imresize, 0=with MTF-filtered downsampling
flag_imresize_MS_RR=0;   % Resizing MS for RR validation 1=with imresize, 0=with MTF-filtered downsampling
flag_imresize_HS_RR=0;   % Resizing HS for RR validation 1=with imresize, 0=with MTF-filtered downsampling
flag_imresize_REF=1;     % Generation of PAN FR reference image 1=with imresize, 0=with MTF-filtered downsampling
flag_imresize_MSREF=0;   % Generation of MS FR reference image 1=with imresize, 0=with MTF-filtered downsampling

%% Methods' parameters
lambda_DMSC=0.001;

%% FR reference image
flag_useREF=0;           % Set image as FR reference image 0=dataset PAN, 1=Multiplatform PAN, 2=Simulated PAN from MS, 3=Multiplatform MS, 4=PAN-ALI, 5= Mean of Multiplatform MS

%% Grouping for HS/MS fusion
flag_grouping=0;         % 0=Groups images according to band overlap; 1=to SAM grouping; 2=to CC grouping; 3=to band distance grouping; 4=to simplified SAM grouping; 5=to SAM grouping with simulated HS; 6= SAM grouping with sim HS (between EXP); 8=CC-group with interp. HS (and Sobel); 9=SAM-based with proper interpolated simulated HS comparison; 10=CC-based with Hist.match, 12=CC-based at FR; 14=SAM-based at FR; 15=SAM-based with HS simulated at FR; 20=estimated centroid distance; 21=estimated dispersion ML; 22=given center distance; 23=given dispersion ML; 24=CCOV-AA; 25=RSR-AA; 26=LSQ-AA (MS); 27=LSQ-AA (HS)
flag_bandbyband=1;       % Calculate band by band SAM, SCC, ERGAS
flag_bandselect_val=0;   % 0=choose all bands for validation; 1=only non overlapping bands; 2=only overlapping bands

%% Quality Indices
flag_RR_val=1;           % 1: classic, 2: my alt, 3: Q2n test, 4: sCC test, 5: FR test at RR, 6: custom
flag_FR_val=5;           % 1: classic, 2: +variants, 3:+Khan & Zhou, 4: +personal, 5: D_s test, 6: Only spatial
flag_R2=0;               % 1: Calculates R2 between PAN and the regression of the low resolution image; 0: doesn't calculate
flag_histmatchREF=2;     % if 1, the reference REF is histogram matched to MS; 2=normalized by radiometric resolution; 3=by max value; 0=non-modified REF
flag_noautoblock=0;      % if 1, forces a value for Q block size
Qblocks_size_setup=32;   % Only valid if above is 1
Qblocks_shift_setup=2;   % Block shift for D_s validation (only used if flag_noautoblock=1 and flag_FR_val=5)
Qblocks_shiftl_setup=16; % Block shift for D_lambda validation (only used if flag_noautoblock=1 and flag_FR_val=5)
q_setup=1;               % q value to use in D_s validation (only if flag_noautoblock=1 and flag_FR_val=5)

%% Methods chosen for custom list

flag_EXP=1;
flag_PCA=1;
flag_IHS=1;
flag_Brovey=1;
flag_BDSD=0;
flag_GS=0;
flag_GSA=1;
flag_GSA_local=0;
flag_PRACS=0;
flag_HPF=0;
flag_SFIM=1;
flag_Indusion=0;
flag_ATWT=1;
flag_AWLP=0;
flag_ATWT_M2=0;
flag_ATWT_M3=0;
flag_ATWT_mod=1;
flag_MTF_GLP=0;
flag_MTF_GLP_HPM_PP=0;
flag_MF_HG=0;
flag_MTF_GLP_HPM=1;
flag_MTF_GLP_CBD=1;
flag_GFPCA=0;
flag_CNMF=1;
flag_BayesNaive=0;
flag_BayesSparse=0;
flag_HySure=1;
flag_G_BDSD=0;
flag_L_BDSD=0;
flag_C_BDSD=0;
flag_BayesNaive_new=0;
flag_BayesSparse_new=0;
flag_BayesML=0;
flag_BayesNLM=0;
flag_DMSC=0;
flag_CASSI=0;