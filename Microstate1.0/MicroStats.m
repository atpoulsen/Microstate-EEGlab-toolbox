%   [OUTEEG] = MicroStats(EEG,epoch)
%   MicroStats calculates microstate statistics.
%   >> OUTEEG = MicroStats( INEEG, 'key1', 'val1', 'key2', 'val2', ....)
%
%  Note - Early untested version.
%
%  Please cite this toolbox as:
%  Poulsen, A. T., Pedroni, A., &  Hansen, L. K. (unpublished manuscript).
%  Microstate EEGlab toolbox: An introductionary guide.
%
%  Inputs:
%  EEG      - EEG-lab EEG structure (channels x samples (x epochs)) with
%             .microstate.fit.bestLabel (created by MicroFit.m).
%
%  Optional input:
%  'epoch'  - timewindow of analysis (vector of timeframes).
%
%  Outputs:
%  OUTEEG.microstate.stats  - Structure of microstate parameters per trial
%   .Gfp        - Global field power
%   .Occurence  - Occurence of a microstate per s
%   .Duration   - Average duration of a microstate
%   .Coverage   - % of time occupied by a microstate
%   .GEV        - Global Explained Variance of microstate
%   .MspatCorr  - Spatial Correlation between template maps and microstates
%   .TP         - transition probabilities
%
%  OUTEEG.microstate.stats.avgs - Structure of microstate parameters mean /
%                                 standard deviation over trials.
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

function [OUTEEG] = MicroStats(EEG,varargin)
%% Error check and initialisation
if nargin < 1
    help MicroPara;
    return;
end;

OUTEEG = EEG;
com = '';

settings = check_settings(varargin, EEG);

%% Write settings to OUTEEG (overwrites any previous microstate info )
OUTEEG.microstate.stats = settings;


%% prepare the data and arrays
% only analyze the timewindow of interest
mslabels = OUTEEG.microstate.fit.bestLabel(:,settings.epoch);
% prepare arrays (only MGFP can have NANs!)
MGFP= nan(size(OUTEEG.data,3),size(OUTEEG.microstate.fit.MStemplate,2));
MDur=zeros(size(OUTEEG.data,3),size(OUTEEG.microstate.fit.MStemplate,2));
MOcc=zeros(size(OUTEEG.data,3),size(OUTEEG.microstate.fit.MStemplate,2));
TCov=zeros(size(OUTEEG.data,3),size(OUTEEG.microstate.fit.MStemplate,2));
GEV = nan(size(OUTEEG.data,3),size(OUTEEG.microstate.fit.MStemplate,2));
MspatCorr = nan(size(OUTEEG.data,3),size(OUTEEG.microstate.fit.MStemplate,2));
GFP = OUTEEG.microstate.fit.GFP(:,settings.epoch);


%% For each MS class...
for ms = 1:size(OUTEEG.microstate.fit.MStemplate,2)
    for t = 1:size(OUTEEG.data,3)
        
        % Mean GFP
        MGFP(t,ms) = nanmean(GFP(t,mslabels(t,:)==ms));
        
        runvalue = []; runs = [];
        [runvalue, runs] = my_RLE(mslabels(t,:));
        
        % Mean Duration
        if isnan(nanmean(runs(runvalue == ms)))
            MDur(t,ms) = 0;
            MOcc(t,ms) = 0;
        else
            MDur(t,ms) =  nanmean(runs(runvalue == ms)) .* (1000 / EEG.srate);
            % Occurence
            MOcc(t,ms) =  length(runs(runvalue == ms))./length(settings.epoch).* OUTEEG.srate;
        end
        % time coverage
        TCov(t,ms) = (MDur(t,ms) .* MOcc(t,ms))./ 1000;
        
        % Average spatial correlation per Microstate
        MspatCorrTMP = EEG.microstate.fit.spatCorr(:,settings.epoch,t);
        MspatCorr(t,ms) = nanmean(MspatCorrTMP(ms,mslabels(t,:)==ms));
        
        % global explained variance
        GEV(t,ms) = (sum(GFP(t,mslabels(t,:)==ms) .* MspatCorrTMP(ms,mslabels(t,:)==ms)).^2)./ (sum(GFP(t,mslabels(t,:)==ms)).^2);
    end
end

%% Transition Probabilities (as with hmmestimate(states,states);)
for t = 1:size(OUTEEG.microstate.fit.order,1)
    states = OUTEEG.microstate.fit.order{t,:};
    states = states(states ~= 0);
    % prepare output matrix
    numStates = size(OUTEEG.microstate.fit.MStemplate,2);
    tr = zeros(numStates);
    seqLen = length(states);
    % count up the transitions from the state path
    for count = 1:seqLen-1
        tr(states(count),states(count+1)) = tr(states(count),states(count+1)) + 1;
    end
    trRowSum = sum(tr,2);
    % if we don't have any values then report zeros instead of NaNs.
    trRowSum(trRowSum == 0) = -inf;
    % normalize to give frequency estimate.
    TP(:,:,t) = tr./repmat(trRowSum,1,numStates);
end


%% Write to EEG structure
%   per trial:
OUTEEG.microstate.stats.Gfp = MGFP;
OUTEEG.microstate.stats.Occurence = MOcc;
OUTEEG.microstate.stats.Duration = MDur;
OUTEEG.microstate.stats.Coverage = TCov;
OUTEEG.microstate.stats.GEV = GEV;
OUTEEG.microstate.stats.MspatCorr = MspatCorr;
OUTEEG.microstate.stats.TP = TP;

if size(OUTEEG.data,3) > 1;
    % mean parameters
    OUTEEG.microstate.stats.avgs.Gfp = nanmean(MGFP,1);
    OUTEEG.microstate.stats.avgs.Occurence = nanmean(MOcc,1);
    OUTEEG.microstate.stats.avgs.Duration = nanmean(MDur,1);
    OUTEEG.microstate.stats.avgs.Coverage = nanmean(TCov,1);
    OUTEEG.microstate.stats.avgs.GEV = nanmean(GEV,1);
    OUTEEG.microstate.stats.avgs.MspatCorr = nanmean(MspatCorr,1);
    
    % standard deviation of parameters
    OUTEEG.microstate.stats.avgs.stdGfp = nanstd(MGFP,1);
    OUTEEG.microstate.stats.avgs.stdOccurence = nanstd(MOcc,1);
    OUTEEG.microstate.stats.avgs.stdDuration = nanstd(MDur,1);
    OUTEEG.microstate.stats.avgs.stdCoverage = nanstd(TCov,1);
    OUTEEG.microstate.stats.avgs.stdGEV = nanstd(GEV,1);
    OUTEEG.microstate.stats.avgs.stdMspatCorr = nanstd(MspatCorr,1);
end
end

% -------------- helper functions -------------- %
function settings = check_settings(vargs, EEG)
%% Check settings
% Checks settings given as optional inputs for MicroStats.
% Undefined inputs is set to default values.
varg_check = { 'epoch'  'integer'    []         1:size(EEG.data,2)};
settings = finputcheck( vargs, varg_check);
if ischar(settings), error(settings); end; % check for error
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
