%CALC_FITMEAS calculates measures of fit.
%
%  Note - Early untested version.
%
%  Calculates measures of fit as defined in [1]. The notation has been
%  adapted to be consistent with the rest of the Microstate EEGlab toolbox.
%  Murray -> toolbox notation: q -> K; r -> k; n_r -> Nk; n -> N; D_r -> Dk.
%
%  [1] - Murray, M. M., Brunet, D., & Michel, C. M. (2008). Topographic
%        ERP analyses: A step-by-step tutorial review. Brain Topography.
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
%  KL      - Krzanowski-Lai criterion.
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
function [KL, W, CV, GEV] = calc_fitmeas(X,A_all,L_all)
%% Initialisation
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
    GFP = var(X);
    % Calculating and summing GEV over all timepoints for the corresponding active microstates
    map_corr = columncorr(X,A(:,L));
    GEV(K_ind) = sum((GFP.*map_corr).^2) / sum(GFP.^2);
    
    
    %% CV
    CV(K_ind) = (sum(sum(X.^2)) - sum(sum(A(:,L).*X).^2)) / (N*(C-1));
    
    
    %% W
    Dk = nan(K,1);
    Nk = nan(K,1);
    for k = 1:K
        cluster = X(:,L==k);
        Nk(k) = size(cluster,2);
        clstrsq = dot(cluster,cluster,1);
        %sum of pair-wise distance between all maps of cluster k
        Dk(k) = sum(sum(bsxfun(@plus,clstrsq',clstrsq)-2*(cluster'*cluster)));
    end
    W(K_ind) = (1./(2*Nk))' * Dk;
    M(K_ind) = W(K_ind).*K^(2/C); %for KL
end


%% KL (excludes first and last microstate segmentations)
if NKs < 3
    % KL requires at least three segmentations
else
    % d(k)=M(k)-M(k+1), excludes last segmentation.
    d = M(1:end-1)-M(2:end);
    % KL_top=d(k-1)-d(k), starts at k+1, i.e. excludes first segmentation
    KL_top = (d(1:end-1) - d(2:end)); 
    KL_bottom = M(1:end-2); % M(k-1)
    KL(2:end-1) = KL_top ./ KL_bottom; %excluding first and last segmentation
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