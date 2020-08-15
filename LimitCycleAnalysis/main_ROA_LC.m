clear all
yalmip('clear')
%close all

%Initialize structs
sys         = struct();
pinit       = struct();
fns         = struct();
numsets     = struct();
filenames   = struct();
plotsets    = struct();

cd ..
addfolders = genpath('LimitCycleAnalysis');
addpath(addfolders);
cd LimitCycleAnalysis


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%%%%%%%% USER INPUT HERE %%%%%%%%%%

plot_only = 0;              % make this 1 to only plot

verimethod  = 'SESSdegdV';    % algorithmic option, others: 'sssdeg2', 'sessdegd', 'eedegd'
system      = 'Goluskin';   % name replacing 'Template' in 'initializeTemplate.m' required!
sys.trafo   = 2;            % choice of MOC, 1 = classic, 2 = center point
sys.Nhp     = 50;           % nr of tau samples, i.e. hyperplanes

pinit.initVscale = 1e-4;   % scaling of initial Gram mat of V, this needs to be tuned once in the very beginning such that feasible for multiplier search! 

numsets.iteration_max   = 50;   % maximum number of iterations in case conv. crit. is not met before
numsets.sdpsetting1     = sdpsettings('solver','mosek','verbose',0,'cachesolvers',1); % add here prefered solver settings for step 1
numsets.sdpsetting2     = sdpsettings('solver','mosek','verbose',1,'cachesolvers',1); % add here prefered solver settings for step 2
numsets.convCrit        = 1e-2; % min relative increase of objective function value, threshold as convergence criteria for terminating the iteration


% Change the numerical settings for the chosen verification method
if strcmp(verimethod,'SSSfixV')
    numsets.degs.V_dU  = 2;     % Do not change.
    numsets.degs.V_dL  = 2;     % Do not change.
    numsets.degs.s1_dU = 6;     % s1 SOS multiplier max degree
    numsets.degs.s1_dL = 2;     % s1 SOS multiplier min degree
    numsets.degs.s2_dU = 2;     % s2 SOS multiplier max degree
    numsets.degs.s2_dL = 0;     % s2 SOS multiplier min degree
    
    verimethod_fn = strcat(verimethod,'_tv',num2str(numsets.tvgamma));     % Do not change.
    
elseif strcmp(verimethod,'SSSdeg2V')
    numsets.degs.V_dU  = 2;     % Do not change.
    numsets.degs.V_dL  = 2;     % Do not change.
    numsets.degs.s1_dU = 6;     % s1 SOS multiplier max degree
    numsets.degs.s1_dL = 2;     % s1 SOS multiplier min degree
    numsets.degs.s2_dU = 2;     % s2 SOS multiplier max degree
    numsets.degs.s2_dL = 0;     % s2 SOS multiplier min degree
   
    verimethod_fn      = strcat(verimethod,'_tv',num2str(numsets.tvgamma)); % Do not change.
    
elseif strcmp(verimethod, 'SESSdegdV')
    numsets.degs.V_dU  = 2;     % Lyapunov function max degree
    numsets.degs.V_dL  = 2;     % Lyapunov function min degree (min. should be 2)
    numsets.degs.s1_dU = 6;     % s1 SOS multiplier max degree
    numsets.degs.s1_dL = 2;     % s1 SOS multiplier min degree
    numsets.degs.s2_dU = 2;     % s2 SOS multiplier max degree
    numsets.degs.s2_dL = 0;     % s2 SOS multiplier min degree
    numsets.degs.s3_dU = 2;     % s3 SOS multiplier max degree
    numsets.degs.s3_dL = 0;     % s3 SOS multiplier min degree
    
    pinit.initBscale    = 0.1*pinit.initVscale;            % try to keep this fixed, reduce if needed (surrogate has to fit into R)
    numsets.ellipsemat  = eye(2);                          % enter here desired fixed surrogate ellipse
    
    numsets.alpha       = 1;                                        % Do not change.
    numsets.ellipsemat  = numsets.ellipsemat/pinit.initBscale;      % Do not change.
    verimethod_fn = verimethod;                                     % Do not change.
    
elseif strcmp(verimethod, 'EEdegdV')
    numsets.degs.V_dU  = 4;     % Lyapunov function max degree
    numsets.degs.V_dL  = 2;     % Lyapunov function min degree (min. should be 2)
    numsets.degs.s1_dU = 6;     % s1 SOS multiplier max degree
    numsets.degs.s1_dL = 2;     % s1 SOS multiplier min degree
    numsets.degs.s2_dU = 2;     % s2 SOS multiplier max degree
    numsets.degs.s2_dL = 0;     % s2 SOS multiplier min degree
    numsets.degs.s3_dU = 2;     % s3 SOS multiplier max degree
    numsets.degs.s3_dL = 0;     % s3 SOS multiplier min degree
    
    numsets.tvball  = 1;        % 0 for same ellipse on each hyperplane, 1 for hyperplane-varying ellipse
    verimethod_fn   = strcat(verimethod,'_tv',num2str(numsets.tvball)); % Do not change.
    
else 
    error('Choose a valid verification method.')
end

%%%%%%%% END: USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

sys.method = verimethod;
sys.system = system;

numsets.clean_thresh    = 1e-6;
numsets.tvgamma         = 0;
numsets.lfac            = 1e-4; 
numsets.epsi            = 1e-5;
numsets.neg_def_fac     = 1e-4;

filenames.initializationFile = strcat('initialize',system);

if exist(filenames.initializationFile,'file')
    filenames.initializationFile = str2func(filenames.initializationFile);
    [sys,fns,pinit] = filenames.initializationFile(sys,fns,pinit,'ROA');
else
    error('System unknown, please create initialization file!')
end

filenames.resultsdirectory      = 'resultsFilesROA/';
filenames.intermedresults       = strcat('intermediate_results_',system,'_',verimethod_fn,'.mat');
filenames.finalresults          = strcat(filenames.resultsdirectory,'final_results_',system,'_',verimethod_fn,'_TF',num2str(sys.trafo),'_V',num2str(numsets.degs.V_dU),'.mat');
filenames.filenamestoplot{1}    = strcat('final_results_',system,'_',verimethod_fn,'_TF',num2str(sys.trafo),'_V',num2str(numsets.degs.V_dU),'.mat'); %or enter each file manually



if plot_only ==0
    % Start verification
        if sys.udim ==0
            fns.periodinitgain = str2func('devalPeriodicLyapunov');
            filenames.periodinitname = strcat('initialLyapunov/LyapunovGains_',system,'_N',num2str(sys.Nhp),'_TF',num2str(sys.trafo),'.mat');    
        else
            fns.periodinitgain = str2func('devalPeriodicRiccati');
            filenames.periodinitname = strcat('initialLyapunov/RiccatiGains_',system,'_N',num2str(sys.Nhp),'_TF',num2str(sys.trafo),'.mat');
        end

        if exist(filenames.intermedresults,'file') && delete_intermediate_file ==1
            delete(filenames.intermedresults)
        end

    veri_filename = str2func(strcat('veriROA_',verimethod));
    veri_filename(sys,fns,filenames, pinit,numsets);
    delete(filenames.intermedresults)
end 


plotROA_LC(filenames,1)
plotROAvolume_LC(filenames,2)


%%% Optional analysis 
%plotHyperplanes_LC(sys)









