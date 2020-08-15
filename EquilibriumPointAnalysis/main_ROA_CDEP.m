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
%%%%%%%% USER INPUT HERE (for binary questions: 1 = yes, 0 = no %%%%%%%%%%

%%% Choose the analyses to perform %%%
comp_roapce     = 1;    % computation of PCE ROA
comp_roastoch   = 1; 	% retrival of stochastic ROA from PCE ROA 
plot_results    = 1;    % choose to plot results


% system settings
system   = 'ShortPeriod';   % name replacing 'Template' in 'initializeTemplate.m' required!
sys.p    = 2;               %PCE truncation order (p in {0,1,2,3..})


% numerical settings for PCE FB design and ROA analysis:

numsets.initVscale = 1e-4;   % initial scaling of Gram mat of V, this needs to be tuned once in the very beginning such that feasible for multiplier search! 

numsets.sdpsetting1          = sdpsettings('solver','mosek','verbose',0, 'beeponproblem', []); % add here prefered solver settings for step 1 (multiplier search)
numsets.sdpsetting2          = sdpsettings('solver','mosek','verbose',0, 'beeponproblem', []); % add here prefered solver settings for step 2 (V optimization)
numsets.sdpsettingK          = sdpsettings('solver','mosek','verbose',0, 'beeponproblem', []); % add here prefered solver settings for  K optimization step
numsets.iteration_max        = 20;  % maximum number of iterations in case conv. crit. is not met before
numsets.iteration_max_preROA = 8;  % maximum number of iterations for initial ROA size before FB design
numsets.convCrit             = 1e-2; % min relative increase of \bar{R}, threshold as convergence criteria for terminating the iteration
numsets.convCrit_preROA      = 2e-1; % min relative increase of R_0, threshold as convergence criteria for terminating the iteration


% numerical settings for stochastic ROA retrieval:

numsetsRE.iteration_max   = 15; % maximum number of iterations in case conv. crit. is not met before
numsetsRE.sdpsetting1     = sdpsettings('solver','mosek','verbose',1, 'beeponproblem', []); %add here prefered solver settings for step 1
numsetsRE.sdpsetting2     = sdpsettings('solver','mosek','verbose',1, 'beeponproblem', []); %add here prefered solver settings for step 2
numsetsRE.convCrit        = 1e-2;


% multiplier and V degrees for PCE FB design and ROA analysis: 

numsets.degs.V_dU  = 2;     % Lyapunov function max degree
numsets.degs.V_dL  = 2;     % Lyapunov function min degree (min. should be 2)
numsets.degs.s1_dU = 4;     % s1 SOS multiplier max degree
numsets.degs.s1_dL = 2;     % s1 SOS multiplier min degree
numsets.degs.s2_dU = 2;     % s2 SOS multiplier max degree
numsets.degs.s2_dL = 0;     % s2 SOS multiplier min degree

% the following are only needed if input constraints included:
numsets.degs.s3_dU = 2;     % s3 SOS multiplier max degree
numsets.degs.s3_dL = 0;     % s3 SOS multiplier min degree
numsets.degs.s4_dU = 2;     % s4 SOS multiplier max degree
numsets.degs.s4_dL = 0;     % s4 SOS multiplier min degree     
numsets.degs.s5_dU = 0;     % s5 SOS multiplier max degree
numsets.degs.s5_dL = 0;     % s5 SOS multiplier min degree
numsets.degs.s6_dU = 0;     % s6 SOS multiplier max degree
numsets.degs.s6_dL = 0;     % s6 SOS multiplier min degree
numsets.degs.s7_dU = 2;     % s7 SOS multiplier max degree
numsets.degs.s7_dL = 0;     % s7 SOS multiplier min degree
numsets.degs.s8_dU = 2;     % s8 SOS multiplier max degree
numsets.degs.s8_dL = 0;     % s8 SOS multiplier min degree
numsets.degs.s9_dU = 0;     % s9 SOS multiplier max degree
numsets.degs.s9_dL = 0;     % s9 SOS multiplier min degree
numsets.degs.s10_dU = 0;    % s10 SOS multiplier max degree
numsets.degs.s10_dL = 0;    % s10 SOS multiplier min degree

numsets.initKvals = 0;  % either insert an initial value in the shape of the actual K, or use a scalar
numsets.initt = 0.1;    % either insert an initial value in the shape of the actual t, or use a scalar
        

% multiplier and V degrees for stochastic ROA retrieval:

numsetsRE.degs.V_dU  = numsets.degs.V_dU;   % Do not change
numsetsRE.degs.V_dL  = numsets.degs.V_dL;   % Do not change
numsetsRE.degs.s1_dU = 0;                   % s1 SOS multiplier max degree
numsetsRE.degs.s1_dL = 0;                   % s1 SOS multiplier min degree
numsetsRE.degs.s2_dU = 2;                   % s2 SOS multiplier max degree
numsetsRE.degs.s2_dL = 0;                   % s2 SOS multiplier min degree
numsetsRE.degs.hi_dU = 2;                   % hi poly multiplier max degree
numsetsRE.degs.hi_dL = 0;                   % hi poly multiplier min degree

numsetsRE.initQ0scale = 1e-3;      % initial scaling of Gram mat of R_0, this needs to be tuned once for initial feasibility of the stochastic ROA retrieval computations (decrease upon infeasibility)

%%%%%%%% END: USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


numsets.gfac                 = 1e-4; %
numsets.clean_thresh         = 1e-6;
numsets.epsi                 = 0;
numsetsRE.clean_thresh       = numsets.clean_thresh;


filenames.initializationFile = strcat('initialize',system);

if exist(filenames.initializationFile,'file')
    filenames.initializationFile = str2func(filenames.initializationFile);
    [sys,fns] = filenames.initializationFile(sys,fns);
else
    error('System unknown, please create initialization file!')
end

if sys.include_input_con ==1
    inputcon_tag = 'wIC_';
elseif sys.include_input_con ==0
    inputcon_tag = 'nIC_';
elseif sys.include_input_con ==-1
    inputcon_tag = 'sOL_';
elseif sys.include_input_con ==-2  
    inputcon_tag = 'nOL_';
end

filenames.resultsdirectory         = 'resultsFilesPCEstoch_CD/';

filenames.intermedresultsROAPCE_pre    = strcat('intermediate_results_PCE_preROA_',system,'_p',num2str(sys.p),'.mat');
filenames.finalresultsROAPCE_pre       = strcat(filenames.resultsdirectory,'final_results_PCE_preROA_',system,'_V',num2str(numsets.degs.V_dU),'_p',num2str(sys.p),'.mat');


filenames.intermedresultsROAPCE = strcat('intermediate_results_PCE_',inputcon_tag,system,'_p',num2str(sys.p),'.mat');
filenames.finalresultsROAPCE     = strcat(filenames.resultsdirectory,'final_results_PCE_',inputcon_tag,system,'_V',num2str(numsets.degs.V_dU),'_p',num2str(sys.p),'.mat');

filenames.finalresultsROAstoch     = strcat(filenames.resultsdirectory,'final_results_',inputcon_tag,system,'_V',num2str(numsets.degs.V_dU),'_var',num2str(sys.varfix(1,1)),'.mat');

filenames.filenamestoplot{1}       = strcat('final_results_'); %or enter each file manually
%filenames.filenamestoplot{1}       = strcat('_V',num2str(numsets.degs.V_dU));
filenames.filenamestoplot{2}       = system;
filenames.filenamestoplot{3}       = 'var0';
% Start verification
if comp_roapce == 1

    veri_filename = str2func('veriROA_PCE_CD');
    veri_filename(sys,fns,filenames,numsets);
        delete(filenames.intermedresultsROAPCE_pre)
    delete(filenames.intermedresultsROAPCE)


end    
if comp_roastoch == 1 
   
    veri_filename = str2func('recover_ROAstoch');
    veri_filename(sys,filenames,numsetsRE);
end


if plot_results ==1
    plotROA_stochEP_CD(filenames)
end


