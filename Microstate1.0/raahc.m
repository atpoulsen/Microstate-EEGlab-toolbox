%RAAHC randomised Atomize and Agglomerate Hierarchical Clustering
%  Performs hierachical clustering interchangingly by either agglomerating
%  or atomizing one cluster at a time. The ratio at which either method is
%  chosen can be set. The default settings are set to the Atomize and
%  Agglormerate Hierarchical Clustering (AAHC) method as described in [1,
%  2].
%  When agglomerating, the two most similar clusters are merged, by
%  measuring inner product (angle) between cluster prototypes. When
%  atomizing, the worst cluster is removed, and each its members assigned
%  to the cluster whos prototype they are most correlated with as described
%  in [1,2,3].
%  Since atomizing using correlation (TAAHC) is sensitive to how it is
%  initialised the order of samples are shuffled to make every start
%  different, making it possible to find the best clustering of multiple
%  runs of the algorithm. 
%
%  To start either AAHC [1] or Topograhpical Atomize and Agglormerate
%  Hierarchical Clustering (TAAHC, as described in [2,3]) use these
%  settings:
%   * AAHC : opts.atom_ratio = 1; opts.atom_measure = 'GEV'.
%   * TAAHC: opts.atom_ratio = 1; opts.atom_measure = 'corr'.
%
%  [1] - Murray, M. M., Brunet, D., & Michel, C. M. (2008). Topographic
%        ERP analyses: A step-by-step tutorial review. Brain Topography.
%  [2] - Poulsen, A. T., Pedroni, A., Langer, N., &  Hansen, L. K.
%        (unpublished manuscript). Microstate EEGlab toolbox: An
%        introductionary guide.
%  [3] - Brunet, D.(2011). Cartool reference Guide. Cartool v3.51.
%
%  Please cite this toolbox as:
%  Poulsen, A. T., Pedroni, A., Langer, N., &  Hansen, L. K. (unpublished
%  manuscript). Microstate EEGlab toolbox: An introductionary guide.
%
%  Andreas Trier Poulsen, atpo@dtu.dk,
%  Lars Kai Hansen.
%  Technical University of Denmark, Cognitive systems - October 2017
%
%  Inputs
%  X       - EEG (channels x samples).
%  K_range - Number of microstates to save prototypes and lables from. Can
%            be a vector.
%  
%  Optional inputs
%  opts. 
%       atom_ratio   - Ratio ([0;1]) between atomizing or agglomerating
%                      clusters. Setting to 1 means exclusively atomizing;
%                      0 exclusively agglomerating (default: 1). 
%       atom_measure - Measure used for finding the worst cluster when
%                      atomizing. 'GEV' uses Global explained variance,
%                      'corr' uses sum squared correlations, and 'varex'
%                      uses variance explained (default: 'corr').
%       verbose      - Print status messages to command window?
%                      (default:1).
%
%  Outputs
%  A_all - Cell containing the spatial distribution of microstate
%          prototypes, for each number of microstates defined in K_range.
%          The dimensions of each cell is (channels x K).
%  L_all - Cell containing the labels of the most active microstate at each
%          timepoint, for each number of microstates defined in K_range.
%          The dimensions of each cell is (1 x samples).

function [A_all, L_all] = raahc(X,K_range,opts)
%% Checking optional inputs
if nargin < 3
   opts = []; 
end
if nargin < 2
   K_range = 1:size(X,2); 
end
if isfield(opts,'atom_ratio')
    atom_ratio = opts.atom_ratio;
else
    atom_ratio = 1;
end
if isfield(opts,'atom_measure') 
    atom_measure = opts.atom_measure;
else
    atom_measure = 'corr';
end
if isfield(opts,'verbose') 
    verbose = opts.verbose;
else
    verbose = 1;
end


%% Initialising
X = double(X); % eigs only support double precision
[C,N] = size(X);
Kmin = min(K_range);

% shuffle X to add randomness to which clusters are removed first (for TAAHC).
shuff_idx = randperm(N);
Xshuff = X(:,shuff_idx);

% preallocate cells for cluster prototypes and labels for all K's to be saved
N_K = length(K_range);
A_all = cell(N_K,1);
L_all = cell(N_K,1);

% Cluster prototype and labels
A = Xshuff./(ones(C,1)*sqrt(sum(Xshuff.*Xshuff,1))); % cluster means (normalised unit length)
K = N;
Lshuff = 1:N;

% initialise R assigment matrix (K x N)
R = speye(N); % sparse assignment matrix

% Calculate GFP for GEV atomise measure
if strcmp(atom_measure, 'GEV')
    GFP2 = std(Xshuff).^2;
    GFP2sum = sum(GFP2);
end


%% Save means and prototypes if K is in K_range
K_ind = find(K == K_range);
if ~isempty(K_ind)
    A_all{K_ind} = A;
    
    % save labels to the corresponding unshuffled X
    L = nan(1,N);
    L(shuff_idx) = Lshuff;
    L_all{K_ind} = L;
end


%% Start removing clusters
if verbose
    fprintf('Starting hierarchical clustering from %i to %i clusters.\n ', K, Kmin)
    tstart = tic;
end
while K >= Kmin    
    if rand>atom_ratio % agglomerate if rand is above atom_ratio
        %% Agglomerate
        if verbose
            fprintf('Agglomerating from %i to %i clusters... ', K, K-1)
        end
        % use inner product (angle) as a polarity invariant similarity measure. 
        cluster_angles = abs(A'*A);
        % set selfsimilarity of clusters to zero (we are interested in similarity between pairs)
        cluster_angles = cluster_angles + eps;
        cluster_angles(eye(K)==1)=0;
        
        % sort angles to find closest clusters
        [angle_sort,sort_idx] = sort(cluster_angles,2,'descend');
        
        % agglomerate two closest clusters 
        [~,k1] = max(angle_sort(:,1));   % first member of closest pair
        k2 = sort_idx(k1,1); % second member of closest pair
        R(k1,:) = R(k1,:) + R(k2,:);   % joint membership of the two clusters
        R = R(1:end~=k2,:); % remove 2nd cluster by contracting the R matrix
        K = size(R,1); % update K
        
        %sanity check that we did not loose members
        Nin=full(sum(sum(R)));
        if N>Nin
            disp(['lost members, ', int2str(Nin)])
            error('Lost member')
        end
        
        % compute new cluster prototypes and labels 
        [Lshuff,~] = find(R); % labels
        A = nan(C,K);
        for k = 1:K
            k_idx = Lshuff==k;
            if sum(k_idx) == 1 
                A(:,k) = Xshuff(:,k_idx)/sqrt(Xshuff(:,k_idx)'*Xshuff(:,k_idx));
            else
                % Use 1st PC => polarity invariant, like Pacual-Marqui (1995).
                Sigk=Xshuff(:,k_idx)*Xshuff(:,k_idx)';
                [Ak, ~] = eigs(Sigk,1);
                A(:,k)=Ak;
            end
        end
        
        if verbose
            fprintf(' Total time spent: %i s.\n ', round(toc(tstart)))
        end
        
    else
        %% Atomize
        if verbose
            fprintf('Atomizing from %i to %i clusters ', K, K-1)
        end
        % compute measure of fit
        switch atom_measure % high values mean good clusters
            case 'GEV'
                if verbose, fprintf('using global explained variance...'), end
                fitmeas = nan(K,1);
                for k = 1:K
                    k_idx = Lshuff==k;
                    fitmeas(k) = sum( GFP2(k_idx)/GFP2sum .* corr(Xshuff(:,k_idx),A(:,k)).^2 );
                end
            case 'corr'
                if verbose, fprintf('using summed squared correlation...'), end
                fitmeas = nan(K,1);
                for k = 1:K
                    k_idx = Lshuff==k;
                    fitmeas(k) = sum(corr(Xshuff(:,k_idx),A(:,k)).^2);
                end
            case 'varex'
                % compute variance explained
                if verbose, fprintf('using variance explained...'), end
                Xrec = zeros(size(Xshuff));
                fitmeas=zeros(K,1);
                for k = 1:K
                    k_idx=full(R(k,:)==1);
                    
                    % find the location of state k
                    UA=A(:,k)*A(:,k)';
                    Xrec(:,k_idx)=UA*Xshuff(:,k_idx);
                    
                    % compute atomize index
                    fitmeas(k)=1 - sum(sum((Xshuff(:,k_idx)-Xrec(:,k_idx)).*(Xshuff(:,k_idx)-Xrec(:,k_idx))))/sum(sum(Xshuff(:,k_idx).*Xshuff(:,k_idx)));
                end
            otherwise
                error('Unknown distance measure selected: %s',atom_measure)
        end
        
        % find members of worst cluster to atomize
        [~,atomC_idx] = min(fitmeas); % index of worst cluster
        atomX_idx = find(R(atomC_idx,:)==1); % indices of members of cluster
        
        % find clusters for lost members (polarity invariant)
        if strcmp(atom_measure,'varex') 
            similarity = abs(A'*Xshuff(:,atomX_idx)); 
        elseif sum(strcmp(atom_measure,{'GEV','corr'}))
            similarity = corr(A, Xshuff(:,atomX_idx)).^2;
        end
        % avoid reassigning to atomized cluster
        similarity = similarity + eps;
        similarity(atomC_idx,:) = 0;  
        
        % reassign lost members to most similar clusters
        [~,sim_idx] = max(similarity,[],1);
        for k=1:K
            R(k,atomX_idx(sim_idx==k)) = 1;
        end
        
        %sanity check that we did not loose members
        Nin=sum(sum(R));
        if N>Nin
            disp(['lost members A', int2str(Nin)])
            error('Lost member')
        end
        
        % removing the worst cluster by contracting the R matrix
        R = R(1:end~=atomC_idx,:);
        
        K = size(R,1); % update K
        
        % check everybody is assigned
        Nin=sum(sum(R));
        if N>Nin
            disp(['lost members B', int2str(Nin)])
        end
        
        % compute new cluster prototypes and labels 
        [Lshuff,~] = find(R); % labels
        A = nan(C,K);
        for k = 1:K
            k_idx = Lshuff==k;
            if sum(k_idx) == 1 
                A(:,k) = Xshuff(:,k_idx)/sqrt(Xshuff(:,k_idx)'*Xshuff(:,k_idx));
            else
                % Use 1st PC => polarity invariant, like Pacual-Marqui (1995).
                Sigk=Xshuff(:,k_idx)*Xshuff(:,k_idx)';
                [Ak, ~] = eigs(Sigk,1);
                A(:,k)=Ak;
            end
        end
        
        if verbose
            fprintf(' Total time spent: %i s.\n ', round(toc(tstart)))
        end
        
    end
    
    %% Save means and prototypes if K is in K_range
    K_ind = find(K == K_range); 
    if ~isempty(K_ind)
        A_all{K_ind} = A;
        
        % save labels to the corresponding unshuffled X
        L = nan(1,N);
        L(shuff_idx) = Lshuff;
        L_all{K_ind} = L;
    end   
end

end

