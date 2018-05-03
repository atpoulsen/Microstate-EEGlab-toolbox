%MICROSMOOTH Temporal smoothing of microstate labels
% Temporal smoothing of microstate labels using either windowed smoothing
% [1] or rejection of smaller segments as described in [2].
% These implementation use the same terminology as used in ICA literature,
% i.e. X = AZ + e, and in [2].
% The implementation of windowed smoothing follows the steps as given in
% [1].
% 
% Usage:
%   >> L = MicroSmooth(X, A, smooth_type,opts)
%
% Please cite this toolbox as:
% Poulsen, A. T., Pedroni, A., Langer, N., &  Hansen, L. K. (2018).
% Microstate EEGlab toolbox: An introductionary guide. bioRxiv.
%
% Inputs:
%  X           - EEG (channels x samples (x trials)).
%  A           - Spatial distribution of microstate prototypes (channels x 
%                K).
%  smooth_type - Smoothing type: 'reject segments' (default) or 'windowed'.
%  opts        - Method-specific settings (struct).
%
% Method-specific inputs:
% * Reject segments:
%   opts.
%        minTime  - Redristibute segments smaller than minTime (in samples)
%                   to the next best fitting microstate (default = 3).
%        polarity - Account for polarity when calculating the global map
%                   dissimilarity. Typically off for spontaneous EEG and on
%                   for ERP data (default = 0).
% * Windowed smoothing:
%   opts.
%        b              - Smoothing width. Integer denoting the number of
%                         samples on each side of current sample
%                         (default: 3).
%        lambda         - Smoothing weight (default: 5).
%        max_iterations - Maximum number of iterations of algorithm
%                         (default: 1000).
%        thresh         - Threshold of convergence based on relative change
%                         in noise variance (default: 1e-6). 
%
% Output:
%  L  - Label of the most active microstate at each timepoint (trials x
%       time).
%
%
%  [1] - Pascual-Marqui, R. D., Michel, C. M., & Lehmann, D. (1995).
%        Segmentation of brain electrical activity into microstates: model
%        estimation and validation. IEEE Transactions on Biomedical
%        Engineering.
%  [2] - Poulsen, A. T., Pedroni, A., Langer, N., &  Hansen, L. K.
%        (unpublished manuscript). Microstate EEGlab toolbox: An
%        introductionary guide.
%
% Authors:
% Andreas Trier Poulsen, atpo@dtu.dk
% Technical University of Denmark, DTU Compute, Cognitive systems.
%
% Andreas Pedroni, andreas.pedroni@uzh.ch
% University of Zürich, Psychologisches Institut, Methoden der
% Plastizitätsforschung.
%
% September 2017.
%
% See also: eeglab pop_micro_segment

% Copyright (C) 2017  Andreas Trier Poulsen, atpo@dtu.dk
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
function L = MicroSmooth(X, A, smooth_type, opts)
%% Initialisation
[C,N,T] = size(X); 
L = nan(T,N);

%% Select smoothing type and loop over trials
switch smooth_type
    case 'windowed'
        for tr = 1:T
            L(tr,:) = window_smoothing(X(:,:,tr),A,opts);
        end
    case 'reject segments'
        for tr = 1:T
            L(tr,:) = reject_segments(X(:,:,tr),A,opts);
        end
    otherwise
    error('Unknown smoothing type: %s', smooth_type)
end
end

% ------------------------- Smoothing functions --------------------------% 
function [L,sig2,ind] = window_smoothing(X,A,opts)
% [L,sig2,R2,MSE,ind] = window_smoothing(X,A,opts)
%  Implementation of the Segmentation Smoothing Algorithm, as described in
%  Table II of [1]. Smoothes using the interval t-b to t+b excluding t.
%  Note, that temporary allocation of labels (denoted with Lambda in [1])
%  is not necessary in this implementation, and steps 3 and 6 are therefore
%  left out.

%% Initialisation (step 1 to 4)
[C,N] = size(X);
K = size(A,2);
const = sum(sum(X.^2));

% Reading settings
if isfield(opts,'lambda')
    lambda = opts.lambda;
else
    lambda = 5;
end
if isfield(opts,'b')
    b = opts.b;
else
    b = 3;
end
if isfield(opts,'max_iterations')
    max_iterations = opts.max_iterations;
else
    max_iterations = 1000;
end
if isfield(opts,'thresh')
    thresh = opts.thresh;
else
    thresh = 1e-6;
end

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
e = (const - sum(sum(A(:,L).*X).^2)) / (N*(C-1));

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
    sig2 = (const - sum(sum(A(:,L).*X).^2)) / (N*(C-1));
    
end

% Step 10
% sig2_D = const / (N*(C-1));
% R2 = 1 - sig2/sig2_D;
% activations = zeros(size(Z));
% for n=1:N; activations(L(n),n) = Z(L(n),n); end % setting to zero
% MSE = mean(mean((X-A*activations).^2));
end

function L = reject_segments(X,A,opts)
% Reject small segments
% If there is a segment of TFs that is smaller than minTime it gets the next
% best label. This starts with segments of length 1 and then iterates up to
% minTime (in samples).

%% Initialisation
[C,N] = size(X);
K = size(A,2);

% Reading settings
if isfield(opts,'minTime')
    minTime = opts.minTime;
else
    minTime = 3;
end
if isfield(opts,'polarity')
    polarity = opts.polarity;
else
    polarity = 0;
end


%% Calculate the Global map dissimilarity for each Microstate
% Normalise EEG and maps (average reference and gfp = 1 for EEG)
X = (X - repmat(mean(X,1), C, 1)) ./ repmat(std(X,1), C, 1);
A = (A - repmat(mean(A,1), C, 1)) ./ repmat(std(A,1), C, 1);

% Global map dissilarity
GMD = nan(K,N);
for k = 1:K
    GMD(k,:) = sqrt(mean( (X - repmat(A(:,k),1,N)).^2 ));
end

% Account for polarity (recommended 0 for spontaneous EEG)
if polarity == 0
    GMDinvpol = nan(K,N);
    for k = 1:K
        GMDinvpol(k,:) = sqrt(mean( (X - repmat(-A(:,k),1,size(X,2))).^2));
    end
    idx = GMDinvpol < GMD;
    GMD(idx) = GMDinvpol(idx);
end

% Sort the GMD to get the labels
[~ ,labels] = sort(GMD,1);


%% reject small maps
for k = 1:minTime
    cruns = k;
    while sum(cruns <= k) > 0
        idx = [];
        [~, runs] = my_RLE(labels(1,:));
        idx(cumsum([1 runs(:,runs>0)])) = 1;
        cruns = runs(:, cumsum(idx(1:find(idx,1,'last')-1)));
        labels(:,cruns<=k) = circshift(labels(:,cruns<=k),-1);
    end
end


%% Export best labels only
L = labels(1,:);
end

function [d,c]=my_RLE(x)
%% RLE
% This function performs Run Length Encoding to a strem of data x.
% [d,c]=rl_enc(x) returns the element values in d and their number of
% apperance in c. All number formats are accepted for the elements of x.
% This function is built by Abdulrahman Ikram Siddiq in Oct-1st-2011 5:15pm.
if nargin~=1
    error('A single 1-D stream must be used as an input')
end

ind=1;
d(ind)=x(1);
c(ind)=1;

for i=2 :length(x)
    if x(i-1)==x(i)
        c(ind)=c(ind)+1;
    else ind=ind+1;
        d(ind)=x(i);
        c(ind)=1;
    end
end
end