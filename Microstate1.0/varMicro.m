function [A_opt,L_opt,Res] = varMicro(X,K_range,opts)
%VARMICRO Variational Microstates.
%
%  Note - Early untested version.
%
%  Implementation of variational microstates, with optional smoothing of
%  segments as described in [1].
%  Free energy is used to choose amongst multiple restarts, for the same 
%  number of microstates (K). Free energy is also used to find the optimal
%  number of microstates.
%
%  [1] - (unpublished manuscript). Variational microstate analysis.
%  [2] - Pascual-Marqui, R. D., Michel, C. M., & Lehmann, D. (1995).
%        Segmentation of brain electrical activity into microstates: model
%        estimation and validation. IEEE Transactions on Biomedical
%        Engineering.
%
%  Please cite this toolbox as:
%  Poulsen, A. T., Pedroni, A., Langer, N., &  Hansen, L. K. (unpublished
%  manuscript). Microstate EEGlab toolbox: An introductionary guide.
%
%  Andreas Trier Poulsen, atpo@dtu.dk
%  Technical University of Denmark, Cognitive systems - February 2017
%
%  Inputs
%  X       - EEG (channels x samples).
%  K_range - A priori number of microstates. Can be a vector.
%  
%  Optional inputs
%  opts. 
%       sig2_0         - Prior variance of activations (default: average
%                        EEG channel variance).
%       p0             - Probability for having same microstate as last 
%                        timepoint (default: 0). Setting to zero turns
%                        smoothing off.
%       reps           - Number of random initialisations of algorithm
%                        (default: 10).
%       max_iterations - Maximum number of iterations of algorithm
%                        (default: 1000).
%       thresh         - Threshold of convergence based on relative change
%                        in free energy (default: 1e-6). 
%       verbose        - Print status messages to command window?
%                        (default:1).
%
%  Outputs
%  A_opt   - Spatial distribution of microstates (channels x K).
%  L_opt   - Label of the most active microstate at each timepoint (1 x samples).
%  Res     - Struct containing additional results:
%     .Z_all    - Cell containing the microstate activations for each number of
%                 microstates defined in K_range. The dimensions of each cell
%                 is (K x samples).
%     .A_all    - Cell containing the spatial distribution for each number of
%                 microstates defined in K_range. The dimensions of each cell
%                 is (channels x K).
%     .L_all    - Cell containing the labels for each number of microstates
%                 defined in K_range. The dimensions of each cell is
%                 (1 x samples).
%     .S_opt         - Sparsity parameters (K x samples) of optimal solution.
%     .beta          - Noise precision for the best solution of each K in K_range.
%     .sig2_z_opt    - Variance of microstate activations (K x samples) of 
%                      optimal solution.
%     .K_act         - Number of microstates for optimal solution.
%     .Free_energy   - Free energy for the best solution of each K in K_range.
%     .sig2_modk     - Noise variance as calculated in [2] for the best 
%                      solution of each K in K_range.
%     .sig2_modk_mcv - Modified predictive residual variance as calculated 
%                      in [2] for the best solution of each K in K_range.
%     .R2            - Explained variance as calculated in [2] for the best
%                      solution of each K in K_range.
%     .MSE           - Mean squared error for the best solution of each K in K_range.
%     .opts          - Settings used to run the algorithm.


%% Checking inputs and setting default settings, where needed.
if nargin<3; opts = []; end
if nargin<2; K_range = 3:8; end

if ~isfield(opts,'reps'), opts.reps = 10; end
if ~isfield(opts,'max_iterations'),	opts.max_iterations = 1000; end
if ~isfield(opts,'verbose'), opts.verbose = 1; end
if ~isfield(opts,'thresh'), opts.thresh = 1e-6; end
if ~isfield(opts,'sig2_0'),	opts.sig2_0 = mean(var(X,0,2)); end
if ~isfield(opts,'p0'), opts.p0 = 0; end

if opts.p0==0 % Checking if smoothing is requested
    opts.smooth = 0;
else
    opts.smooth = 1;
end   


%% Preallocating
N_K = length(K_range);
A_all = cell(N_K,1);
L_all = cell(N_K,1);
Z_all = cell(N_K,1);
R2_all = nan(N_K,1);
sig2_modk_all = nan(N_K,1);
sig2_modk_mcv_all = nan(N_K,1);
MSE_all = nan(N_K,1);
beta_all = nan(N_K,1);
F_all = nan(N_K,1);


%% Looping over all K values in K_range
F_opt = inf;
K_ind = 0;
if opts.verbose
    if opts.smooth
        disp('Starting variational microstates with temporal smoothing')
    else
        disp('Starting variational microstates without temporal smoothing')
    end
end

for K = K_range
    K_ind = K_ind + 1;
    if opts.verbose
        fprintf('Analysis no. %i out of %i. Starting %i random initialisations for %i microstates.\n'...
            ,K_ind,N_K,opts.reps, K)
    end
    
    % Finding best fit amongst a given number of restarts
    F_best = inf;
    for r = 1:opts.reps
        if opts.verbose
            fprintf('Starting initialisations no. %i out of %i. '...
                ,r,opts.reps)
        end
        
        % Segmenting with variational microstates
        [A,Z,sig2_z,S,beta,F,ind,sig2_modk,sig2_modk_mcv,R2,MSE] = segmentation(X,K,opts);
        if opts.verbose
            fprintf('Finished in %i iterations.\n',ind)
        end
        
        % Checking for best fit
        if F < F_best
            A_all{K_ind} = A;
            Z_all{K_ind} = Z;
            % Calculating label of active microstate for each timepoint for
            % best solution
            [~,L_all{K_ind}] = max(S);
            sig2_z_best = sig2_z;
            S_best = S;
            beta_all(K_ind) = beta;
            F_best = F;
            sig2_modk_all(K_ind) = sig2_modk;
            sig2_modk_mcv_all(K_ind) = sig2_modk_mcv;
            R2_all(K_ind) = R2;
            MSE_all(K_ind) = MSE;
        end
    end
    
    F_all(K_ind) = F_best;
    
    % Finding optimum solution amongst different values of K
    if F_all(K_ind) < F_opt
        A_opt = A_all{K_ind};
        L_opt = L_all{K_ind};
        sig2_z_opt = sig2_z_best;
        S_opt = S_best;
        K_act = K;
        F_opt = F_all(K_ind);
    end
end


%% Saving additional results and info to Res struct
Res.Z_all = Z_all;
Res.A_all = A_all;
Res.L_all = L_all;
Res.S_opt = S_opt;
Res.beta = beta_all;
Res.sig2_z_opt = sig2_z_opt;
Res.K_act = K_act;
Res.Free_energy = F_all;
Res.sig2_modk = sig2_modk_all;
Res.sig2_modk_mcv = sig2_modk_mcv_all;
Res.R2 = R2_all;
Res.MSE = MSE_all;
Res.opts = opts;

end

function [A,Z,sig2_Z,S,beta,F,ind,sig2_modk,sig2_modk_mcv,R2,MSE] = segmentation(X,K,opts)
% Updating variables for variational microstates in iterations.
% Dimensions of variables: X=(C x N), Z=(K x N), sig2_Z=(K x N), A=(C x K),
% S=(K x N).

%% Initialising
[C,N] = size(X);
F = inf;
F_old = 0;

% Reading settings
max_iterations = opts.max_iterations;
thresh = opts.thresh;
sig2_0 = opts.sig2_0;
p0 = opts.p0;

tmp = (1-p0)/(K-1);
smooth_const = log(1 + (p0 - tmp)/tmp);

% random initialisation of A, beta, and S
A = X(:,randperm(N,K)); % selecting K random timepoints to use as initial microstate maps
A = bsxfun(@rdivide,A,sqrt(sum(A.*A)));% normalising

beta = 1/mean(var(X,0,2));

S = rand(K,N)*0.1; % giving all states a non-zero probability
[~,L] = max((A'*X).^2); % finding most active microstates based on A
for n = 1:N, S(L(n),n) = 1; end % setting these to 1
S = bsxfun(@rdivide,S,sum(S,1));% normalising so each n sum to 1


%% Iterate over updates until convergence
ind = 0;
while abs(F_old-F) >= thresh*F && max_iterations>ind
    ind = ind + 1;
    F_old = F;
    % Z (mu_z)
    denom = 1/(sig2_0) + beta*bsxfun(@times, S,sum(A.^2)');
    Z = beta*(S.*(A'*X)) ./ denom;
    if isnan(Z), error('hov'), end
    
    % sig2_z
    sig2_Z = 1./denom;
    if isnan(sig2_Z), error('hov'), end
    
    % S (mu_s) with optional temporal smoothing
    if opts.smooth
        S_past = [zeros(K,1) S(:,1:end-1)];
        S_future = [S(:,2:end) zeros(K,1)];
    end
    
    S = beta*( Z.*(A'*X) - 0.5*bsxfun(@times, Z.^2+sig2_Z, sum(A.^2)') );
    if isnan(S), error('hov'), end
    if opts.smooth
        S = S + smooth_const*(S_past+S_future);
    end
    S = safe_softmax(S);
    if isnan(S), error('hov'), end
    
    % A
    A = bsxfun(@rdivide, X*(Z.*S)', sum((Z.^2+sig2_Z).*S,2)');
    if isnan(A), error('hov'), end
    
    % beta
    beta = C*N / sum(sum( X.^2 - 2*X.*(A*(Z.*S)) + A.^2*((Z.^2+sig2_Z).*S) ));
    if isnan(beta), error('hov'), end
    
    % Free energy
    F = free_energy(Z,sig2_Z,sig2_0,S,beta,C,N,K);
end
%% Calculating measures of fit
% Noise variance and performance measures as calculated in Pascual-Marqui (1995)
const = sum(sum(X.^2));
[~,L] = max(Z.^2); %finding most active microstates
sig2_modk = (const - sum(sum(A(:,L).*X).^2)) / (N*(C-1));
sig2_modk_mcv = sig2_modk * ((C-1)^-1 * (C-1-K))^-2;
sig2_D = const / (N*(C-1));
R2 = 1 - sig2_modk/sig2_D;

% Mean squared error
MSE = mean(mean((X-A*(Z.*S)).^2));

end

function X = safe_softmax(X)
% Computes a safe softmax, to avoid troubles with numerical precision, by
% subtracting largest values before exponentiation using exp(x) =
% exp(x-a)*exp(a).
% Input
% X (K,N) log-probs
% Output
% X (K,N) softmax'ed values

% Subtracting max values for each n
A = max(X,[],1);
X = bsxfun(@minus,X,A);

% Softmaxing
X = exp(X);
X = bsxfun(@times,X,1./sum(X,1));

% Avoid NaNs when calculating log(X)
X(X==0) = eps; 
end

function F = free_energy(Z,sig2,sig2_0,S,beta,C,N,K)
% Calculate free energy of solution. beta needs to be last variable
% calculated.
F1 = -0.5*(N*K*log(2*pi)+sum(sum(log(sig2)))) -N*K/2 + sum(sum(S.*log(S)));
F2 = 0.5*K*N*log(2*pi*sig2_0) + sum(sum(Z.^2+sig2))/(2*sig2_0) + N*log(K);
F3 = -C*N/2*log(beta/(2*pi)) + C*N/2; % using expression for beta in last part

F = F1+F2+F3;
if isnan(F), error('hov'), end
end