%CALC_FITMEAS calculates measures of fit.
%
%  Calculates measures of fit as defined in [1]. The notation has been
%  adapted to be consistent with the rest of the Microstate EEGlab toolbox.
%  Murray -> toolbox notation: q -> K; r -> k; n_r -> Nk; n -> N; D_r -> Dk.
%  Includes the Krzanowski-Lai criterion as originally described in [2] 
%  with the added rule that KL(K)=0 if the dispersion function increases.
%  
%
%  [1] - Murray, M. M., Brunet, D., & Michel, C. M. (2008). Topographic
%        ERP analyses: A step-by-step tutorial review. Brain Topography.
%  [2] - Krzanowski, W., & Lai, Y. (1988). A criterion for determining the
%        number of groups in a dataset using sum of squares clustering.
%        Biometrics.
%
%  Please cite this toolbox as:
%  Poulsen, A. T., Pedroni, A., Langer, N., &  Hansen, L. K. (unpublished
%  manuscript). Microstate EEGlab toolbox: An introductionary guide.
%
%  Inputs:
%  X       - EEG (channels x samples).
%  A_all   - Cell array of spatial distributions (channels x K), for
%            different number of microstates. Can also be single 2D array.
%  L_all   - Cell array of labels of the most active microstate at each
%            timepoint (1 x samples). Can also be single 1D array.
%
%  Outputs:
%  KL      - Krzanowski-Lai criterion as described in [2].
%  KL_nrm  - Normalised Krzanowski-Lai criterion as described in [1].
%  W       - Dispersion of clusters.
%  CV      - Cross-validation criterion.
%  GEV     - Global explained variance.
%
%  Authors:
%  Andreas Trier Poulsen, atpo@dtu.dk
%  Franciszek Zdyb, follzd@dtu.dk
%  Technical University of Denmark, DTU Compute, Cognitive systems.
%
%  February 2017.

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
function [KL, KL_nrm, W, CV, GEV] = calc_fitmeas(X,A_all,L_all)
%% Initialisation
GFP = std(X);

if iscell(A_all)
    NKs = length(A_all);
else
    NKs = 1;
end

[C,N] = size(X);
GEV = nan(NKs,1);
CV = nan(NKs,1);
W = nan(NKs,1);
M = nan(NKs,1); 
KL = nan(NKs,1);
KL_nrm = nan(NKs,1);


%% Looping over K's
for K_ind = 1:NKs
    %% Obtaining mpas and labels for K microstates
    if iscell(A_all)
        A = A_all{K_ind};
        L = L_all{K_ind};
    else
        A = A_all;
        L = L_all;
    end
    K = size(A,2);
    
    
    %% GEV
    % Calculating and summing GEV over all timepoints for the corresponding active microstates
    map_corr = columncorr(X,A(:,L));
    GEV(K_ind) = sum((GFP.*map_corr).^2) / sum(GFP.^2);
    
    
    %% CV
    nrm = ( (C-1)/(C-1-K) )^2;
    CV(K_ind) = nrm * (sum(sum(X.^2)) - sum(sum(A(:,L).*X).^2)) / (N*(C-1));
    
    
    %% W
    Dk = nan(K,1);
    Nk = nan(K,1);
    for k = 1:K
        cluster = X(:,L==k);
        Nk(k) = size(cluster,2);
        clstrsq = dot(cluster,cluster,1);
        %sum of pair-wise distance between all maps of cluster k
        Dk(k) = sum(sum(bsxfun(@plus,clstrsq',clstrsq)-2*(cluster'*cluster)));
        
        % polarity invariant
%         cluster_corr = abs(corr(cluster));
%         cluster_corr_pairs = triu(cluster_corr,1);
%         Dk2(k) = mean(cluster_corr_pairs(:));
        
    end
    idx = Nk ~= 0; % in case of empty clusters
    W(K_ind) = (1./(2*Nk(idx)))' * Dk(idx);
    M(K_ind) = W(K_ind)*K^(2/C); %for KL
    
    % polarity invariant
%     W2(K_ind) = mean(Dk2);
%     M2(K_ind) = W2(K_ind)*K^(2/C); %for KL
    
end


%% KL_nrm - Normalised KL as described in [1]
% (excludes first and last microstate segmentations)
if NKs < 3
    % KL requires at least three segmentations
else
    % preallocating
    d = nan(NKs,1);
    KL_top = nan(NKs,1);
    KL_bottom = nan(NKs,1);
    
    % d(K)=M(K)-M(K+1), excludes last segmentation.
    d(1:end-1) = M(1:end-1) - M(2:end);
    
    % KL_top=d(K-1)-d(K), starts at K+1, i.e. excludes first segmentation
    KL_top(2:end) = d(1:end-1) - d(2:end);
    KL_bottom(2:end) = M(1:end-1);

    KL_nrm = KL_top ./ KL_bottom; 
    
    % only convex shapes of the W is considered for KL. Note in [1] they
    % mistankenly write concave.
    % KL(K) = 0 if d(K-1)<0 or d(K-1)<d(K)
    idx1 = [false; d(1:end-1)<0];
    idx2 = [false; d(1:end-1)<d(2:end)];
    KL_nrm(idx1) = 0;
    KL_nrm(idx2) = 0;
    KL_nrm([1 end]) = nan;
end


%% KL - KL as described in [2]
% (excludes first and last microstate segmentations)
if NKs < 3
    % KL requires at least three segmentations
else
    % preallocating
    diff = nan(NKs,1);
    
    % diff(K)=M(K-1)-M(K), excludes first segmentation. note: different from KL_nrm.
    diff(2:end) = M(1:end-1) - M(2:end);
    
    % KL=abs(diff(K)/diff(K+1)), excludes last segmentation
    KL(1:end-1) = abs(diff(1:end-1) ./ diff(2:end));
    
    % Added rule that W(K) - W(K-1) cannot be positive, i.e. W increases from K-1
    % to K.    
    idx_KL = [false; W(2:end) - W(1:end-1)];
    KL(idx_KL>0) = 0;
    KL([1 end]) = nan;
end


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