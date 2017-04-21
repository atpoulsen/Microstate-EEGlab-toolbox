%  [OUTEEG] = MicroFit(EEG,MStemp,varargin)
%
%  MicroFit fits microstate prototype maps to data.
%  >> OUTEEG = MicroFit (EEG, MStemp, 'key1', 'val1', 'key2', 'val2', ....)
%
%  Note - Early untested version.
%
%  Please cite this toolbox as:
%  Poulsen, A. T., Pedroni, A., &  Hansen, L. K. (unpublished manuscript).
%  Microstate EEGlab toolbox: An introductionary guide.
% 
%  Inputs
%  EEG              - EEG-lab EEG structure with EEG.data (channels x time (x trials)).
%  MStemp           - Microstates Template maps (channels x N microstates)
%
% Optional inputs:
%  'minTF'          - Redristibute segments smaller than minTF to the next best
%                     fitting microstate (default = 0)
%  'polarity'       - account for polarity when fitting (spontaneous EEG
%                     typically ignore polarity = 0, ERP data = 1) (default = 0)
%  'sequentialize'  - sequentialize clusters (spontaneous EEG
%                     typically ignore sequence microstates = 0, ERP data = 1) (default = 0)
%  'getorder'       - get the order of microstates required for MicroPara
%                     computation of Transition Probabilities (default = 1)
%
%  Outputs
%  OUTEEG.microstate.fit.mslabels       - Microstate labels (N microstates x time (x trials))
%  OUTEEG.microstate.fit.gmd            - GMD between N Microstate Prototype maps (N  x time (x trials))
%  OUTEEG.microstate.fit.seq            - Sequence of reoccurrence of MS ((trials x)  time)
%  OUTEEG.microstate.fit.bestLabel      - best fitting MS per TF ((trials x)  time)
%  OUTEEG.microstate.fit.spatCorr       - spatial correlation Template / EEG (N  x time (x trials))
%  OUTEEG.microstate.fit.msFirstTF      - first occurence of a microstate (similar to use like seq) ((trials x)  time)
%  OUTEEG.microstate.fit.minTF          - see inputs
%  OUTEEG.microstate.fit.polarity       - see inputs
%  OUTEEG.microstate.fit.sequentialize  - see inputs
%  OUTEEG.microstate.fit.MStemplate     - see inputs
%
% Author:
%
% Andreas Pedroni, andreas.pedroni@uzh.ch
% University of Zürich, Psychologisches Institut, Methoden der
% Plastizitätsforschung.
%
% February 2017.

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

function [OUTEEG] = MicroFit(EEG,MStemp,varargin)

%% Error check and initialisation
if nargin < 2
    help MicroFit;
    return;
end;

OUTEEG = EEG;

% force average reference
OUTEEG.data = OUTEEG.data - repmat(mean(OUTEEG.data,1),[size(OUTEEG.chanlocs,2),1,1]);

settings = check_settings(varargin);

%% Write settings to OUTEEG (overwrites any previous microstate.fit info )
OUTEEG.microstate.fit = settings;

%% check if the data has more than one trial or not and reshape if necessary
if size(EEG.data,3) > 1;
    X = squeeze(reshape(EEG.data,size(EEG.data,1),size(EEG.data,2)*size(EEG.data,3)));
else
    X = EEG.data;
end

%% calculate the Global map dissimilarity for each Microstates

% a average reference and gfp = 1 for both maps
X = (X - repmat(mean(X,1),size(EEG.data,1),1)) ./ repmat(std(X,1),size(EEG.data,1),1);
MStemp = (MStemp - repmat(mean(MStemp,1),size(EEG.data,1),1)) ./ repmat(std(MStemp,1),size(EEG.data,1),1);

% b dissimilarity
for k = 1:size(MStemp,2)
    GMD1(k,:) = sqrt(mean( (X - repmat(MStemp(:,k),1,size(X,2))).^2));
    GMD2(k,:) = sqrt(mean( (X - repmat(-MStemp(:,k),1,size(X,2))).^2));
end


% if spontaneous EEG do not account for polarity of maps
if settings.polarity == 0
    take1 = GMD2 >= GMD1;
    take2 = GMD2 < GMD1;
    GMD = take1 .* GMD1 + take2 .* GMD2;
else
    GMD = GMD1;
end

% Sort the GMD to get the labels
[~ ,label] = sort(GMD,1);

%% reject small segments
% if there is a segment of TFs that is smaller than minTF it gets the next
% best label. This starts with segments of length 1 and then iterates up to
%  minTF
for k = 1:settings.minTF
    cruns = k;
    while sum(cruns <= k) > 0
        idx = [];
        [~, runs] = my_RLE(label(1,:));
        idx(cumsum([1 runs(:,runs>0)])) = 1;
        cruns = runs(:, cumsum(idx(1:find(idx,1,'last')-1)));
        label(:,cruns<=k) = circshift(label(:,cruns<=k),-1);
        % do the same for GMD...
        GMD(:,cruns<=k) = circshift(GMD(:,cruns<=k),-1);
    end
end


%% calculate the spatial correlation between Template map and EEG
SpatCorr = 1 - (GMD.^2)./2;
SpatCorr = squeeze(reshape(SpatCorr,size(MStemp,2),size(EEG.data,2),size(EEG.data,3)));

%% get the Labels of the best fitting Microstates.
bestLabel = label(1,:);


%% Sequentialize
% take into account the sequence of occurence of microstates. This makes
% only sense in ERP data (if it makes sense)
if settings.sequentialize == 1
    if size(EEG.data,3) > 1
        bestLabel = squeeze(reshape(bestLabel,1,size(EEG.data,2),size(EEG.data,3)))';
    end
    
    seq = zeros(size(bestLabel,1),size(bestLabel,2));
    for t = 1:size(bestLabel,1)
        first = 1;
        s = ones(size(MStemp,2),1);
        for k = 1:size(bestLabel,2)-1
            if bestLabel(t,k) == bestLabel(t,k+1)
                seq(t,k) = s(bestLabel(t,k))  ;
                msFirstTF(t,k) = first;
            else
                seq(t,k) = s(bestLabel(t,k))  ;
                s(bestLabel(t,k),1) = s(bestLabel(t,k),1) + 1;
                first = k;
                msFirstTF(t,k) = first;
            end
        end
        seq(t,k+1) = seq(t,k);
        msFirstTF(t,k+1) = msFirstTF(t,k);
    end
else
    % this needs to be done for the transition probabilities
    if size(EEG.data,3)  > 1
        bestLabel = reshape(bestLabel,size(EEG.data,2),size(EEG.data,3))';
    end
end
% % quick plot to check
% subplot(3,1,1)
% area(bestLabel(1,:))
% subplot(3,1,2)
% area(seq(1,:))
% subplot(3,1,3)
% area(msFirstTF(1,:))




%% MS Order (for transition probabilities in MicroPara)
if settings.getorder == 1
    order = {};
    for t = 1:size(bestLabel,1)
        [order{t,1}, ~ ] = my_RLE(bestLabel(t,:));
    end
end


%% add GFP to output
GFP = squeeze(std(OUTEEG.data,[],1))';

% correct the number of trials in EEG
%% only for ERP!
% add Prototype maps
if size(EEG.data,3)  == 1
    OUTEEG.microstate.fit.mslabels = label;
    OUTEEG.microstate.fit.bestLabel = bestLabel;
    OUTEEG.microstate.fit.gmd = GMD;
    OUTEEG.microstate.fit.spatCorr = SpatCorr;
    OUTEEG.microstate.fit.GFP = GFP';
    if settings.sequentialize == 1
        OUTEEG.microstate.fit.seq = seq;
        OUTEEG.microstate.fit.msFirstTF = msFirstTF;
    end
    if settings.getorder == 1;
        OUTEEG.microstate.fit.order = order;
    end
    
else
    OUTEEG.microstate.fit.mslabels = reshape(label,size(label,1),size(EEG.data,2),size(EEG.data,3));
    OUTEEG.microstate.fit.bestLabel = bestLabel ;
    OUTEEG.microstate.fit.gmd = squeeze(reshape(GMD,size(GMD,1),size(EEG.data,2),size(EEG.data,3)));
    OUTEEG.microstate.fit.spatCorr = SpatCorr;
    OUTEEG.microstate.fit.GFP = GFP;
    if settings.sequentialize == 1
        OUTEEG.microstate.fit.seq = seq;
        OUTEEG.microstate.fit.msFirstTF = msFirstTF;
    end
    if settings.getorder == 1;
        OUTEEG.microstate.fit.order = order;
    end
end
% add Microstate Template map
OUTEEG.microstate.fit.MStemplate = MStemp;
end

% -------------- helper functions -------------- %

function settings = check_settings(vargs)
%% check settings
% Checks settings given as optional inputs for MicroFit.
% Undefined inputs is set to default values.
varg_check = {  'minTF'         'integer'	[]	0;
    'polarity'      'integer'	[]	0;
    'sequentialize' 'integer'	[]	0;
    'getorder'      'integer'	[]	1};

settings = finputcheck( vargs, varg_check);
if ischar(settings), error(settings); end; % check for error
end

function [d,c]=my_RLE(x)
%% RLE
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
