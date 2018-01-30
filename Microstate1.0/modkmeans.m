function [A_opt,L_opt,Res] = modkmeans(X,K_range,opts)
%MODKMEANS Modified K-means algorithm.
%  Implementation of the Modified K-means algorithm, with optional
%  smoothing of segments as described in [1]. This implementation use the
%  same terminology as used in ICA litetrature, i.e. X = AZ + e. The
%  implementation follows the steps as given in [1]. Note that the
%  non-active microstates in "Z" are not set zero.
%  A new and optimised method for segmentation is also available. This
%  method and upholds the same model assumptions, but has a different 
%  iteration scheme. Preliminary tests show that this method is faster and
%  is even slightly better at representing EEG with microstates. See [2]
%  for more information.
%  Explained variance is used to choose amongst multiple restarts, after
%  smoothing, for the same number of microstates (K). To find the
%  optimal number of microstates the modified predictive residual variance,
%  sig2_mcv from [1] is used.
%
%  [1] - Pascual-Marqui, R. D., Michel, C. M., & Lehmann, D. (1995).
%        Segmentation of brain electrical activity into microstates: model
%        estimation and validation. IEEE Transactions on Biomedical
%        Engineering.
%  [2] - Poulsen, A. T., Pedroni, A., Langer, N., &  Hansen, L. K.
%        (unpublished manuscript). Microstate EEGlab toolbox: An
%        introductionary guide.
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
%       b              - Smoothing width. Denotes the samples on each side
%                        of current sample (default: 0). Setting to zero
%                        turns smoothing off.
%       lambda         - Smoothing weight (default: 5).
%       reps           - Number of random initialisations of algorithm
%                        (default: 10).
%       max_iterations - Maximum number of iterations of algorithm
%                        (default: 1000).
%       thresh         - Threshold of convergence based on relative change
%                        in noise variance (default: 1e-6).
%       verbose        - Print status messages to command window?
%                        (default:1).
%       fitmeas        - Readying measure of fit for selecting best
%                        segmentation. 'CV': Crossvalidation criterion,
%                        'GEV': Global explained variance, 'dispersion':
%                        Dispersion of clusters. (Default: 'CV').
%       optimised      - Use the new and optimised segmentation introduced
%                        in [2]? (Default: 0).
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
%     .R2       - Explained variance for the best solution of each K in K_range.
%     .sig2_modk     - Noise variance for the best solution of each K in K_range.
%     .sig2_modk_mcv - Modified predictive residual variance of best solution
%                      of each K in K_range.
%     .MSE      - Mean squared error for the best solution of each K in K_range.
%     .K_act    - Number of microstates for optimal solution.
%     .opts     - Settings used to run the algorithm.


%% Initialising and checking optional inputs
const1 = sum(sum(X.^2));
C = size(X,1);

% Checking inputs and setting default settings, where needed.
if nargin<3; opts = []; end
if nargin<2; K_range = 3:8; end

if ~isfield(opts,'reps'), opts.reps = 10; end
if ~isfield(opts,'max_iterations'),	opts.max_iterations = 1000; end
if ~isfield(opts,'verbose'), opts.verbose = 1; end
if ~isfield(opts,'thresh'), opts.thresh = 1e-6; end
if ~isfield(opts,'lambda'),	opts.lambda = 5; end
if ~isfield(opts,'b'), opts.b = 0; end
if ~isfield(opts,'fitmeas'), opts.fitmeas = 'CV'; end
fitmeas = opts.fitmeas;
if ~isfield(opts,'optimised'), opts.optimised = 0; end

if opts.b==0 % Checking if smoothing is requested
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
MSE_all = nan(N_K,1);
sig2_all = nan(N_K,1);


%% Readying measure of fit for selecting best segmentation
switch fitmeas
    case 'GEV'
        GFP = std(X);
        GFP_const = sum(GFP.^2);
        GEV_opt = 0;
    case 'CV'
        sig2_mcv_opt = inf;
    case 'dispersion'
        W_opt = inf;
    otherwise
        error('Unknown measuref of fit requested: %s',fitmeas)
end
new_best = 0;


%% Looping over all K values in K_range
K_ind = 0;

if opts.verbose
    if opts.optimised
        algo_str = 'the new and optimised ';
    else
        algo_str = '';
    end
    if opts.smooth
        smooth_str = 'with';
    else
        smooth_str = 'without';
    end
    fprintf('Starting %smodified K-means %s temporal smoothing.\n',...
        algo_str, smooth_str)
end

for K = K_range
    K_ind = K_ind + 1;
    if opts.verbose
        fprintf('Analysis no. %i out of %i. Starting %i random initialisations for %i microstates.\n'...
            ,K_ind,N_K,opts.reps, K)
    end
    
    % Finding best fit amongst a given number of restarts based on selected
    % measure of fit
    switch fitmeas
        case 'GEV'
            GEV_best = 0;
        case 'CV'
            sig2_all(K_ind) = inf;
        case 'dispersion'
            W_best = inf;
    end
    
    for r = 1:opts.reps
        if opts.verbose
            fprintf('Starting initialisations no. %i out of %i. '...
                ,r,opts.reps)
        end
        
        
        if opts.optimised
            % The new and optimised Basic N-Microstate Algorithm (see [2])
            [A,L,Z,sig2,R2,MSE,ind] = opt_segment(X,K,const1,opts);
        else
            % The original Basic N-Microstate Algorithm (Table I in [1])
            [A,L,Z,sig2,R2,MSE,ind] = segmentation(X,K,const1,opts);
        end
        
        
        if opts.verbose
            fprintf('Finished in %i iterations. ',ind)
        end
        
        % The Segmentation Smoothing Algorithm (Table II)
        if opts.smooth
            if opts.verbose
                fprintf('Smoothing... ')
            end
            [L,sig2,R2,MSE,ind] = smoothing(X,A,K,const1,opts);
            if opts.verbose
                fprintf('Smoothing done in %i iterations.\n',ind)
            end
        else
            if opts.verbose, fprintf('\n'), end
        end
        
        % Checking for best fit.
        switch fitmeas
            case 'GEV'
                map_corr = columncorr(X,A(:,L));
                GEV = sum((GFP.*map_corr).^2) / GFP_const;
                if GEV > GEV_best, new_best = 1; GEV_best = GEV; end
            case 'CV'
                % We use sig2 instead of R2 as they are proportional
                if sig2 < sig2_all(K_ind), new_best = 1; sig2_all(K_ind) = sig2; end
            case 'dispersion'
                W = calc_dispersion(X,L,K);
                if W < W_best, new_best = 1; W_best = W; end
        end
        
        if new_best
            A_all{K_ind} = A;
            L_all{K_ind} = L;
            sig2_all(K_ind) = sig2;
            Z_all{K_ind} = Z;
            R2_all(K_ind) = R2;
            MSE_all(K_ind) = MSE;
            new_best = 0;
        end
    end
    
    sig2_mcv(K_ind) = sig2_all(K_ind) * ((C-1)^-1 * (C-1-K))^-2;
    
    % Checking for best fit for different values of K.
    switch fitmeas
        case 'GEV'
            if GEV_best > GEV_opt, new_best = 1; GEV_opt = GEV_best; end
        case 'CV'
            % We use sig2 instead of R2 as they are proportional
            if sig2_mcv(K_ind) < sig2_mcv_opt, new_best = 1; sig2_mcv_opt = sig2_mcv(K_ind); end
        case 'dispersion'
            if W_best < W_opt, new_best = 1; W_opt = W_best; end
    end
    
    if new_best
        A_opt = A_all{K_ind};
        L_opt = L_all{K_ind};
        K_act = K;
        new_best = 0;
    end
    
end


%% Saving additional results and info to Res struct
Res.Z_all = Z_all;
Res.A_all = A_all;
Res.L_all = L_all;
Res.R2 = R2_all;
Res.MSE = MSE_all;
Res.sig2_modk = sig2_all;
Res.sig2_modk_mcv = sig2_mcv;
Res.K_act = K_act;
Res.opts = opts;

end

function [A,L,Z,sig2,R2,MSE,ind] = segmentation(X,K,const1,opts)
%  Implementation of the Basic N-Microstate Algorithm, as described in
%  Table I of [1].

%% Initialising (step 1 and 2a 3)
[C,N] = size(X);

% Reading settings
max_iterations = opts.max_iterations;
thresh = opts.thresh;

% Step 1
sig2_old = 0;
sig2 = Inf;

% Step 2a
A = X(:,randperm(N,K)); % selecting K random timepoints to use as initial microstate maps
A = bsxfun(@rdivide,A,sqrt(diag(A*A')));% normalising

%% Iterations (step 3 to 6)
ind = 0;
while abs(sig2_old-sig2) >= thresh*sig2 && max_iterations>ind
    ind = ind + 1;
    sig2_old = sig2;
    
    % Step 3
    Z = A'*X;
    [~,L] = max(Z.^2);
    
    % Step 4
    for k = 1:K
        k_idx = L==k;
        if sum(k_idx)==0
            % no members of this microstate
            A(:,k) = 0;
        else
            S = X(:,k_idx)*X(:,k_idx)';
            %finding eigenvector with largest value and normalising it
            [eVecs,eVals] = eig(S,'vector');
            [~,idx] = max(abs(eVals));
            A(:,k) = eVecs(:,idx);
            A(:,k) = A(:,k)./sqrt(sum(A(:,k).^2));
        end
    end
    
    % Step 5
    sig2 = (const1 - sum(sum(A(:,L).*X).^2)) / (N*(C-1));
    
end


%% Saving solution converged on (step 7 and 8)
% Step 7
Z = A'*X; % NOTE, not setting non-activated microstates to zero
[~,L] = max(Z.^2);

% Step 8
sig2_D = const1 / (N*(C-1));
R2 = 1 - sig2/sig2_D;
activations = zeros(size(Z));
for n=1:N; activations(L(n),n) = Z(L(n),n); end % setting to zero
MSE = mean(mean((X-A*activations).^2));
end

function [A,L,Z,sig2,R2,MSE,ind] = opt_segment(X,K,const1,opts)
%  New and optimised implementation of the Basic N-Microstate Algorithm,
%  as described in [2].

%% Initialising (step 1 and 2a 3)
[C,N] = size(X);
deltaZ = nan(N,1); % preallocating

% Reading settings
max_iterations = opts.max_iterations;
thresh = opts.thresh;

% Step 1
sig2_old = 0;
sig2 = Inf;

% Step 2a
A = X(:,randperm(N,K)); % selecting K random timepoints to use as initial microstate maps
A = bsxfun(@rdivide,A,sqrt(diag(A*A')));% normalising

%% Iterations (step 3 to 6)
ind = 0;
while abs(sig2_old-sig2) >= thresh*sig2 && max_iterations>ind
    ind = ind + 1;
    sig2_old = sig2;
    
    % Find activations
    Z = A'*X;
    [~,L] = max(Z.^2);
    for n=1:N
        deltaZ(n) = Z(L(n),n);
    end 
    
    % Find A
    for k = 1:K
        k_idx = L==k;
        if sum(k_idx)==0
            % no members of this microstate
            A(:,k) = 0;
        else
            A(:,k) = X(:,k_idx)*deltaZ(k_idx);
            A(:,k) = A(:,k)./sqrt(sum(A(:,k).^2));
        end
    end
    
    % Estimate residual noise
    sig2 = (const1 - sum(sum(A(:,L).*X).^2)) / (N*(C-1));
    
end


%% Saving solution converged on (step 7 and 8)
% Step 7
Z = A'*X; % NOTE, not setting non-activated microstates to zero
[~,L] = max(Z.^2);

% Step 8
sig2_D = const1 / (N*(C-1));
R2 = 1 - sig2/sig2_D;
activations = zeros(size(Z));
for n=1:N; activations(L(n),n) = Z(L(n),n); end % setting to zero
MSE = mean(mean((X-A*activations).^2));
end

function [L,sig2,R2,MSE,ind] = smoothing(X,A,K,const1,opts)
%  Implementation of the Segmentation Smoothing Algorithm, as described in
%  Table II of [1]. Smoothes using the interval t-b to t+b excluding t.
%  Note, that temporary allocation of labels (denoted with Lambda in [1])
%  is not necessary in this implementation, and steps 3 and 6 are therefore
%  left out.

%% Initialisation (step 1 to 4)
[C,N] = size(X);

% Reading settings
lambda = opts.lambda;
b = opts.b;
max_iterations = opts.max_iterations;
thresh = opts.thresh;

% Step 1
sig2_old = 0;
sig2 = Inf;

% Step 2
Z = A'*X;
[~,L] = max(Z.^2);

%Check to avoid the loop getting caught and switching one label back and
% forth between iterations.
L_old{1} = zeros(size(L));
L_old{2} = zeros(size(L));

% Step 4
e = (const1 - sum(sum(A(:,L).*X).^2)) / (N*(C-1));

% Defining constant for step 5b
tmp = sum(X.^2);
const_5b = (repmat(tmp,K,1) - Z.^2) / (2*e*(C-1));
clear tmp


%% Iterations (step 5 to 8)
ind = 0;
while abs(sig2_old-sig2) >= thresh*sig2 && max_iterations>ind ...
        && mean(L_old{rem(ind,2)+1} == L)~=1
    ind = ind + 1;
    sig2_old = sig2;
    L_old{abs(rem(ind,2)-2)} = L;
    
    % Step 5a
    Nbkt_tmp = zeros(K,N);
    for k = 1:K
        Nbkt_tmp(k,:) = double(L==k);
    end
    %using filter to count the number of labels equal to k before (tmp1)
    %and after (tmp2) a given timepoint.
    tmp1 = filter([0 ones(1,b)],1,Nbkt_tmp,[],2);
    tmp2 = filter([0 ones(1,b)],1,Nbkt_tmp(:,end:-1:1),[],2);
    Nbkt = tmp1 + tmp2(:,end:-1:1);
    
    % Step 5b
    [~,L] = min( const_5b - lambda*Nbkt );
    
    % Step 7
    sig2 = (const1 - sum(sum(A(:,L).*X).^2)) / (N*(C-1));
    
end

% Step 10
sig2_D = const1 / (N*(C-1));
R2 = 1 - sig2/sig2_D;
activations = zeros(size(Z));
for n=1:N; activations(L(n),n) = Z(L(n),n); end % setting to zero
MSE = mean(mean((X-A*activations).^2));
end

function C2 = columncorr(A,B)
% Fast way to compute correlation of multiple pairs of vectors without
% computing all pairs as would with corr(A,B). Borrowed from Oli at Stack
% overflow. Note the resulting coefficients vary slightly from the ones
% obtained from corr due differences in the order of the calculations.
% (Differences are of a magnitude of 1e-9 to 1e-17 depending of the tested
% data).

An=bsxfun(@minus,A,mean(A,1));
Bn=bsxfun(@minus,B,mean(B,1));
An=bsxfun(@times,An,1./sqrt(sum(An.^2,1)));
Bn=bsxfun(@times,Bn,1./sqrt(sum(Bn.^2,1)));
C2=sum(An.*Bn,1);

end

function W = calc_dispersion(X,L,K)
% Calculate dispersion. Same code as in calc_fitmeas.m
Dk = nan(K,1);
Nk = nan(K,1);
for k = 1:K
    cluster = X(:,L==k);
    Nk(k) = size(cluster,2);
    clstrsq = dot(cluster,cluster,1);
    %sum of pair-wise distance between all maps of cluster k
    Dk(k) = sum(sum(bsxfun(@plus,clstrsq',clstrsq)-2*(cluster'*cluster)));
end
idx = Nk ~= 0; % in case of empty clusters
W = (1./(2*Nk(idx)))' * Dk(idx);
end
