clear all
clear yalmip


sys         = struct();
fns         = struct();
numsets     = struct();
numsetsRE     = struct();
filenames   = struct();

cd ..
addfolders = genpath('EquilibriumPointAnalysis');
addpath(addfolders);
cd EquilibriumPointAnalysis


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%%%%%%%% BEGIN: USER INPUT (for binary questions: 1 = yes, 0 = no) %%%%%%%%


%%% Choose the analyses to perform %%%
comp_roapce     = 1;    % computation of PCE ROA
comp_roastoch   = 1; 	% retrival of stochastic ROA from PCE ROA 
plot_results    = 1;    % choose to plot results

% system settings
system      = 'VDPEP';  % name replacing 'Template' in 'initializeTemplate.m' required!
sys.p       = 3;        % PCE truncation order (p in {0,1,2,3..})

% numerical settings for PCE ROA analysis:
numsets.initVscale = 1e-4;   % initial scaling of Gram mat of V, this needs to be tuned once in the beginning such that feasible for multiplier search (decrease upon infeasibility)

numsets.iteration_max   = 50;   % maximum number of iterations in case conv. crit. is not met before
numsets.sdpsetting1     = sdpsettings('solver','mosek','verbose',1); %add here prefered solver settings for step 1
numsets.sdpsetting2     = sdpsettings('solver','mosek','verbose',1); %add here prefered solver settings for step 2
numsets.convCrit        = 1e-2; % min relative increase of \bar{R}, threshold as convergence criteria for terminating the iteration


% numerical settings for stochastic ROA retrieval:

numsetsRE.iteration_max   = 15; % maximum number of iterations in case conv. crit. is not met before
numsetsRE.sdpsetting1     = sdpsettings('solver','mosek','verbose',1); %add here prefered solver settings for step 1
numsetsRE.sdpsetting2     = sdpsettings('solver','mosek','verbose',1); %add here prefered solver settings for step 2
numsetsRE.convCrit        = 1e-2; % min relative increase of R_0, threshold as convergence criteria for terminating the iteration


delete_intermediate_file_roapce = 0; % turn on only for initial tuning


% multiplier and V degrees for PCE ROA analysis: 

numsets.degs.V_dU  = 4;            % Lyapunov function max degree
numsets.degs.V_dL  = 2;            % Lyapunov function min degree (min. should be 2)
numsets.degs.s1_dU = 4;            % s1 SOS multiplier max degree 
numsets.degs.s1_dL = 2;            % s1 SOS multiplier min degree
numsets.degs.s2_dU = 2;            % s2 SOS multiplier max degree
numsets.degs.s2_dL = 0;            % s2 SOS multiplier min degree


% multiplier and V degrees for stochastic ROA retrieval:

numsetsRE.degs.V_dU  = numsets.degs.V_dU;   % Do not change
numsetsRE.degs.V_dL  = numsets.degs.V_dL;   % Do not change
numsetsRE.degs.s1_dU = 0;                   % s1 SOS multiplier max degree  
numsetsRE.degs.s1_dL = 0;                   % s1 SOS multiplier min degree
numsetsRE.degs.s2_dU = 2;                   % s2 SOS multiplier max degree
numsetsRE.degs.s2_dL = 0;                   % s2 SOS multiplier min degree
numsetsRE.degs.hi_dU = 2;                   % hi poly multiplier max degree
numsetsRE.degs.hi_dL = 0;                   % hi poly multiplier min degree

numsetsRE.initQ0scale = 1e-3;    % initial scaling of Gram mat of R_0, this needs to be tuned once for initial feasibility of the stochastic ROA retrieval computations (decrease upon infeasibility)

filenames.filenamestoplot{1}      = strcat('final_results_',system,'_V',num2str(numsets.degs.V_dU)); %adjust to filename of desired results to plot

%%%%%%%% END: USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

numsets.gfac            = 1e-4; 
numsets.clean_thresh    = 1e-6; 
numsetsRE.clean_thresh  = numsets.clean_thresh;

filenames.initializationFile = strcat('initialize',system);

if exist(filenames.initializationFile,'file')
    filenames.initializationFile = str2func(filenames.initializationFile);
    [sys,fns] = filenames.initializationFile(sys,fns);
else
    error('System unknown, please create initialization file!')
end

filenames.resultsdirectory        = 'resultsFilesPCEstoch/';
filenames.intermedresultsROAPCE   = strcat('intermediate_results_PCE_',system,'_p',num2str(sys.p),'.mat');
filenames.finalresultsROAPCE      = strcat(filenames.resultsdirectory,'final_results_PCE_',system,'_V',num2str(numsets.degs.V_dU),'_p',num2str(sys.p),'.mat');

filenames.finalresultsROAstoch    = strcat(filenames.resultsdirectory,'final_results_',system,'_V',num2str(numsets.degs.V_dU),'_var',num2str(sys.varfix(1,1)),'.mat');


if comp_roapce == 1

    veri_filename = str2func('veriROA_PCE_');
    veri_filename(sys,fns,filenames,numsets);

    save(filenames.finalresultsROAPCE, '-append','numsets','sys')
    delete(filenames.intermedresultsROAPCE)
end    

if comp_roastoch == 1

    veri_filename = str2func('recover_ROAstoch');
    veri_filename(sys,filenames,numsetsRE);

    save(filenames.finalresultsROAstoch, '-append','numsetsRE','sys')
end


if plot_results ==1
    plotROA_stochEP(filenames)
end









