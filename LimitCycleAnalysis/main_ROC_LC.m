clear all
yalmip('clear')

%Initialize structs
sys             = struct();
pinit           = struct();
fns             = struct();
roc.numsets     = struct();
uncer.numsets   = struct();
filenames       = struct();

cd ..
addfolders = genpath('LimitCycleAnalysis');
addpath(addfolders);
cd LimitCycleAnalysis


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%%%%%%%% USER INPUT HERE %%%%%%%%%%

plot_only = 0;                  % make this 1 to only plot

verimethod  = 'maxUncertaintyInv';      % analysis type: 'maxROCinv', 'maxUncertaintyInv', (use maxROCballInv if invariant region should not be metric sublevel set)
system      = 'QuadPlanLC';     % name replacing 'Template' in 'initializeTemplate.m' required!
sys.trafo   = 1;                % choice of MOC, 1 = classic, 2 = center point
sys.Nhp     = 50;               % nr of tau samples, i.e. hyperplanes

% ------
% Note: a 'fixed ROC estimate' for the maxUncertainty analysis is here
% obtained by performing a rough roa maximization in the beginning. If data
% of a desired shape is available rename the file to the name in filenames.rocfinalresults
% ------


% fill these out for both choices, 'maxROC' and 'maxUncertainty'

roc.numsets.iteration_max   = 50; % make this small if 'maxUncertainty' chosen, otherwise allow for enough iterations for precise convergence
roc.numsets.convCrit        = 1e-2; % make this small if 'maxROC' chosen, otherwise this can be larger to go quicker to the uncertainty max
roc.numsets.sdpsetting1     = sdpsettings('solver','mosek','verbose',0,'cachesolvers',1); %add here prefered solver settings for step 1
roc.numsets.sdpsetting2     = sdpsettings('solver','mosek','verbose',0,'cachesolvers',1); %add here prefered solver settings for step 2

roc.numsets.degs.M_dU   = 2;    % Metric max degree
roc.numsets.degs.M_dL   = 0;    % Metric min degree
roc.numsets.degs.sm1_dU = 4; 	% m1 matrix SOS multiplier max degree
roc.numsets.degs.sm1_dL = 0;    % m1 matrix SOS multiplier min degree
roc.numsets.degs.s2_dU  = 2;    % s2 SOS multiplier max degree 
roc.numsets.degs.s2_dL  = 0;    % s2 SOS multiplier min degree 
roc.numsets.degs.s3_dU  = 2;    % s3 SOS multiplier max degree  
roc.numsets.degs.s3_dL  = 0;    % s3 SOS multiplier min degree 
roc.numsets.degs.h1_dU  = 4;    % h1 indef multiplier max degree 
roc.numsets.degs.h1_dL  = 0;    % h1 indef multiplier min degree 

pinit.initVscale    = 5e-5;                      % scaling of initial metric, this needs to be tuned once in the very beginning such that feasible for multiplier search! 
pinit.initBscale    = 0.1*pinit.initVscale;      % try to keep this fixed, reduce if needed (surrogate has to fit into Z)
ellipsemat          = eye(2);                    % enter here desired fixed surrogate ellipse

roc.numsets.ellipsemat  = ellipsemat/pinit.initBscale;      % Do not change.
roc.numsets.tvalpha     = 1;                                % 0 for same sublevel set size on each hyperplane, 1 for hyperplane-varying sizes
roc.numsets.alpha       = 1;                                % Do not change. 


% fill these out if 'maxUncertainty' is chosen

uncer.numsets.iteration_max     = 50; % max number of iterations in case convergence criteria is not met before
uncer.numsets.sdpsetting1       = sdpsettings('solver','mosek','verbose',0,'cachesolvers',1); %add here prefered solver settings for step 1
uncer.numsets.sdpsetting2       = sdpsettings('solver','mosek','verbose',0,'cachesolvers',1); %add here prefered solver settings for step 2
uncer.numsets.sdpsetting3       = sdpsettings('solver','mosek','verbose',0,'cachesolvers',1); %add here prefered solver settings for step 3
uncer.numsets.convCrit          = 5e-3; % min relative increase of uncertainty bound, threshold as convergence criteria for terminating the iteration

uncer.numsets.degs = roc.numsets.degs;  % recommended, otherwise copy paste the multipliers above, replace 'roc.' by 'uncer.' and fill in desired values

uncer.numsets.metric_search = 1;    % keep metric fixed or allow the search for more suitable ones


%%%%%%%% END: USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


filenames.resultsdirectory = 'resultsFilesROC/';

filenames.initializationFile = strcat('initialize',system);

if exist(filenames.initializationFile,'file')
    filenames.initializationFile = str2func(filenames.initializationFile);
    [sys,fns,pinit] = filenames.initializationFile(sys,fns,pinit,'ROC');
else
    error('System unknown, please create initialization file!')
end


if (strcmp(verimethod,'maxROCinv') || strcmp(verimethod,'maxROCballInv'))
    
    verimethod_fn = strcat(verimethod,'_tv',num2str(roc.numsets.tvalpha));
    
    filenames.rocintermedresults   = strcat('intermediate_results_',system,'_',verimethod_fn,'_TF',num2str(sys.trafo),'.mat');
    filenames.rocfinalresults      = strcat(filenames.resultsdirectory,'final_results_',system,'_',verimethod_fn,'_TF',num2str(sys.trafo),'_M',num2str(roc.numsets.degs.M_dU),'_var',num2str(sys.uncermax),'.mat');
   
    uncer.numsets = [];

elseif strcmp(verimethod,'maxUncertaintyInv') 
    
    filenames.rocintermedresults   = strcat('intermediate_results_',system,'_initial',verimethod,'_tv',num2str(roc.numsets.tvalpha),'_TF',num2str(sys.trafo),'.mat');
    filenames.rocfinalresults      = strcat(filenames.resultsdirectory,'final_results_',system,'_initial',verimethod,'_tv',num2str(roc.numsets.tvalpha),'_TF',num2str(sys.trafo),'_M',num2str(roc.numsets.degs.M_dU),'.mat');
       
    verimethod_fn = strcat(verimethod,'_tv',num2str(roc.numsets.tvalpha),'_Mvar',num2str(uncer.numsets.metric_search));
    
    filenames.uncerintermedresults   = strcat('intermediate_results_',system,'_',verimethod_fn,'_TF',num2str(sys.trafo),'.mat');
    filenames.uncerfinalresults      = strcat('resultsFilesROC/final_results_',system,'_',verimethod_fn,'_TF',num2str(sys.trafo),'_M',num2str(roc.numsets.degs.M_dU),'.mat');
    
else 
    error('Choose a valid verification method.')
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

sys.method = verimethod;
sys.system = system;

roc.numsets.clean_thresh    = 1e-6;
uncer.numsets.clean_thresh  = roc.numsets.clean_thresh;     


% Start verification
if sys.udim ==0
    fns.periodinitgain = str2func('devalPeriodicLyapunov');
    filenames.periodinitname = strcat('initialLyapunov/LyapunovGains_',system,'_N',num2str(sys.Nhp),'_TF',num2str(sys.trafo),'.mat');    
else
    fns.periodinitgain = str2func('devalPeriodicRiccati');
    filenames.periodinitname = strcat('initialLyapunov/RiccatiGains_',system,'_N',num2str(sys.Nhp),'_TF',num2str(sys.trafo),'.mat');
end

if plot_only == 0 
    if exist(filenames.rocintermedresults,'file') && delete_intermediate_file_maxROC ==1
        delete(filenames.rocintermedresults)
    end
    if strcmp(verimethod, 'maxUncertainty')
        if exist(filenames.uncerintermedresults,'file') && delete_intermediate_file_maxUncer ==1
            delete(filenames.uncerintermedresults)
        end
    end

    veri_filename = str2func(strcat('veriROC_',verimethod));
    veri_filename(sys,fns,filenames, pinit,roc.numsets, uncer.numsets);
end

%delete(filenames.intermedresults)

%filenames.filenamestoplot{1}   = strcat('final_results_',system,'_',verimethod); %or enter each file manually
%filenames.filenamestoplot{2}   = 'TF2';


%plotROC_LC(filenames,10)











