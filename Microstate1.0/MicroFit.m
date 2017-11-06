%MICROFIT Backfitting microstate prototype maps to data.
%
% Usage:
%  >> L = MicroFit(X,A)
%  >> L = MicroFit(X,A,polarity)
%
%  Please cite this toolbox as:
%  Poulsen, A. T., Pedroni, A., &  Hansen, L. K. (unpublished manuscript).
%  Microstate EEGlab toolbox: An introductionary guide.
% 
%  Inputs:
%   X - EEG (channels x samples (x trials)).
%   A - Spatial distribution of microstate prototypes (channels x K).
%
%  Optional input:
%   polarity - Account for polarity when fitting. Typically off for
%              spontaneous EEG and on for ERP data (default = 0).
%
%  Output:
%   L  - Label of the most active microstate at each timepoint (trials x
%        time).
%
% Authors:
%
% Andreas Pedroni, andreas.pedroni@uzh.ch
% University of Zürich, Psychologisches Institut, Methoden der
% Plastizitätsforschung.
%
% Andreas Trier Poulsen, atpo@dtu.dk
% Technical University of Denmark, DTU Compute, Cognitive systems.
%
% September 2017.

% Copyright (C) 2017 Andreas Pedroni, andreas.pedroni@uzh.ch.
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

function L = MicroFit(X,A,polarity)
%% Error check and initialisation
if nargin < 2
    help MicroFit;
    return;
end
if ~exist('polarity','var')
   polarity = 0; 
end

[C,N,T] = size(X);
K = size(A,2);

% force average reference
X = X - repmat(mean(X,1),[C,1,1]);


%% Check if the data has more than one trial or not and reshape if necessary
if T > 1 % for epoched data
    X = squeeze(reshape(X, C, N*T));
end

%% Calculate the Global map dissimilarity for each Microstates
% Normalise EEG and maps (average reference and gfp = 1 for EEG)
X = X ./ repmat(std(X,1), C, 1); % already have average reference
A = (A - repmat(mean(A,1), C, 1)) ./ repmat(std(A,1), C, 1);

% Global map dissilarity
GMD = nan(K,N*T);
for k = 1:K
    GMD(k,:) = sqrt(mean( (X - repmat(A(:,k),1,N*T)).^2 ));
end

% Account for polarity? (recommended 0 for spontaneous EEG)
if polarity == 0
    % Polarity invariance
    GMDinvpol = nan(K,N*T);
    for k = 1:K
        GMDinvpol(k,:) = sqrt(mean( (X - repmat(-A(:,k),1,size(X,2))).^2));
    end
    idx = GMDinvpol < GMD;
    GMD(idx) = GMDinvpol(idx);
end

% Sort the GMD to get the labels
[~ ,labels] = sort(GMD,1);


%% Export best labels only
L = labels(1,:);

if T > 1 % for epoched data
    L = squeeze(reshape(L, 1, N, T))';
end
end