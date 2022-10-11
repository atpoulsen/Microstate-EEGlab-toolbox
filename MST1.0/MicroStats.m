%MICROSTATS Calculates microstate statistics.
%
% Usage:
%   >> Mstats_inclRaw = MicroStats(X, A, L)
%   >> Mstats_inclRaw = MicroStats(X, A, L,polarity,fs)
%
%  Please cite this toolbox as:
%  Poulsen, A. T., Pedroni, A., Langer, N., &  Hansen, L. K. (2018).
%  Microstate EEGlab toolbox: An introductionary guide. bioRxiv.
%
%  Inputs:
%   X - EEG (channels x samples (x trials)).
%   A - Spatial distribution of microstate prototypes (channels x K).
%   L - Label of the most active microstate at each timepoint (trials x
%       time).
%
%  Optional input:
%   polarity - Account for polarity when fitting Typically off for
%              spontaneous EEG and on for ERP data (default = 0).
%   fs       - Sampling frequency of EEG (default = 1).
%
%  Outputs:
%  Mstats - Structure of microstate parameters per trial:
%   .Gfp        - Global field power
%   .Occurence  - Occurence of a microstate per s
%   .Duration   - Average duration of a microstate
%   .Coverage   - % of time occupied by a microstate
%   .GEV        - Global Explained Variance of microstate
%   .MspatCorr  - Spatial Correlation between template maps and microstates
%   .TP         - transition probabilities
%   .seq        - Sequence of reoccurrence of MS ((trials x)  time).
%   .msFirstTF  - First occurence of a microstate (similar to use like seq)
%                 ((trials x)  time)
%   .polarity   - see inputs.
%
%  Mstats.avgs - Structure of microstate parameters mean / and std over
%                trials.
% MTruninger added this 17.02.2021
%  Mstats.raw - Structure of raw microstate parameters* per trial 
%               (for duration, GFP, GEV and spatial correlation)
%               *raw meaning calculated separately for each single
%               occurring microstate
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
%
% Some additions by Moritz Truninger, moritz.truninger@uzh.ch.
%
% February 2021

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

function Mstats = MicroStats(X,A,L,polarity,fs)
%% Error check and initialisation
if nargin < 5
    fs = 1;
elseif nargin < 3
    polarity = 0;
elseif nargin < 3
    help MicroStats;
    return;
end

[C,N,Ntrials] = size(X);
K = size(A,2);
GFP = squeeze(std(X));

%% Spatial correlation between microstate prototypes and EEG
% force average reference
X = X - repmat(mean(X,1),[C,1,1]);

%Check if the data has more than one trial or not and reshape if necessary
if Ntrials > 1 % for epoched data
    X = squeeze(reshape(X, C, N*Ntrials));
    L_temp = L';
    L_vec = reshape(L_temp, 1, N*Ntrials);
end

% Normalise EEG and maps (average reference and gfp = 1 for EEG)
Xnrm = X ./ repmat(std(X,1), C, 1); % already have average reference
A_nrm = (A - repmat(mean(A,1), C, 1)) ./ repmat(std(A,1), C, 1);

% Global map dissilarity
if Ntrials>1 % changed by NL 3.3.2021 GMD = nan(K,N);
GMD = nan(K,N*Ntrials); 
for k = 1:K
    GMD(k,:) = sqrt(mean( (Xnrm - repmat(A_nrm(:,k),1,N*Ntrials)).^2 )); % changed by NL 3.3.2021  GMD(k,:) = sqrt(mean( (Xnrm - repmat(A_nrm(:,k),1,N)).^2 ));
end
else
    GMD = nan(K,N); 
for k = 1:K
    GMD(k,:) = sqrt(mean( (Xnrm - repmat(A_nrm(:,k),1,N)).^2 ));
end
end
    

% Account for polarity (recommended 0 for spontaneous EEG)
if polarity == 0
    GMDinvpol = nan(K,N*Ntrials);
    for k = 1:K
        GMDinvpol(k,:) = sqrt(mean( (Xnrm - repmat(-A_nrm(:,k),1,size(Xnrm,2))).^2));
    end
    idx = GMDinvpol < GMD;
    GMD(idx) = GMDinvpol(idx);
end

% Calculate the spatial correlation between microstate prototypes and EEG
SpatCorr = 1 - (GMD.^2)./2;

% APedroni added this 13.12.2017
% GEVtotal = sum((SpatCorr(sub2ind(size(SpatCorr),L,1:size(L,2))).*squeeze(std(X)).^2) ./sum(squeeze(std(X)).^2));
% MTruninger changed this 15.12.2020
if Ntrials>1 
    GEVtotal = sum(((SpatCorr(sub2ind(size(SpatCorr),L_vec,1:size(L_vec,2))).*squeeze(std(X))).^2) ./sum(squeeze(std(X)).^2));
else
GEVtotal = sum(((SpatCorr(sub2ind(size(SpatCorr),L,1:size(L,2))).*squeeze(std(X))).^2) ./sum(squeeze(std(X)).^2));
end

SpatCorr = squeeze(reshape(SpatCorr,K,N,Ntrials));

%% Sequentialize
% take into account the sequence of occurence of microstates. This makes
% only sense in ERP data (if it makes sense)
seq = nan(Ntrials, N);
msFirstTF = nan(Ntrials, N);
for trial = Ntrials
    first = 1;
    s = ones(K,1);
    for n = 1:(N-1)
        if L(trial,n) == L(trial,n+1)
            seq(trial,n) = s(L(trial,n));
            msFirstTF(trial,n) = first;
        else
            seq(trial,n) = s(L(trial,n))  ;
            s(L(trial,n),1) = s(L(trial,n),1) + 1;
            first = n;
            msFirstTF(trial,n) = first;
        end
    end
    seq(trial,n+1) = seq(trial,n);
    msFirstTF(trial,n+1) = msFirstTF(trial,n);
end


%% Microstate order for transition probabilities
order = cell(Ntrials,1);
for trial = 1:Ntrials
    [order{trial}, ~ ] = my_RLE(L(trial,:));
end


%% Preallocating arrays for stats and readying GFP
% prepare arrays (only MGFP can have NANs!)
MGFP = nan(Ntrials,K);
MDur = zeros(Ntrials,K);
MOcc = zeros(Ntrials,K);
TCov = zeros(Ntrials,K);
GEV = nan(Ntrials,K);
MspatCorr = nan(Ntrials,K);

if Ntrials>1
    GFP = GFP';
end

%% Preallocating structure and arrays for raw stats
% MTruninger added this section 17.02.2021
raw = struct;
for trial2 = 1:Ntrials
    
    % name of current trial (i.e. epoch)
    n_trial = join(["trial",string(trial2)],'');
    
    % get runvalue (= label) of each single microstate and runs
    % (= duration of a single microstate in timepoints)
    [runvalue, runs] = my_RLE(L(trial2,:));
    
    % save microstate sequence (sequence of microstate labels)
    raw.(n_trial).sequence = runvalue;
    
    % get timepoints belonging to each single microstate
    timepoints = tp_indivMS(runs);
    raw.(n_trial).timepoints = timepoints;
    
    % prepare arrays within the structure
    raw.(n_trial).GFP = zeros(1, length(runvalue));
    raw.(n_trial).Duration = zeros(1, length(runvalue));
    raw.(n_trial).MspatCorr = zeros(1, length(runvalue));
    raw.(n_trial).GEV = zeros(1, length(runvalue));
    
    clear runvalue runs
end

%% For each MS class...
for k = 1:K
    for trial = 1:Ntrials
        
        % Mean GFP
        MGFP(trial,k) = nanmean(GFP(trial,L(trial,:)==k));
        
        [runvalue, runs] = my_RLE(L(trial,:));
        
        % Mean Duration
        if isnan(nanmean(runs(runvalue == k)))
            MDur(trial,k) = 0;
            MOcc(trial,k) = 0;
        else
            MDur(trial,k) =  nanmean(runs(runvalue == k)) .* (1000 / fs);
            % Occurence
            MOcc(trial,k) =  length(runs(runvalue == k))./N.* fs;
        end
        % time coverage
        TCov(trial,k) = (MDur(trial,k) .* MOcc(trial,k))./ 1000;
        
        % Average spatial correlation per Microstate
        MspatCorrTMP = SpatCorr(:,:,trial);
        MspatCorr(trial,k) = nanmean(MspatCorrTMP(k,L(trial,:)==k));
        
        % global explained variance Changed by Pedroni 3.1.2018
        GEV(trial,k) = sum( (GFP(trial,L(trial,:)==k) .* MspatCorrTMP(k,L(trial,:)==k)).^2) ./ sum(GFP(trial,:).^2);
        
    end
end

%% Same for raw stats...
% MTruninger added this section 17.02.2021
% (this could possibly be implemented in the part above)
for k = 1:K %for each MS class...
    for trial = 1:Ntrials
     
        [runvalue, runs] = my_RLE(L(trial,:));
            
        % name of current trial (i.e. epoch)
        n_trial = join(["trial",string(trial)],'');
        
        % needed below
        MspatCorrTMP = SpatCorr(:,:,trial);
        
        for tt = find(raw.(n_trial).sequence(:)'==k)
           
            % Raw GFP
            raw.(n_trial).GFP(tt) = ...
                nanmean(GFP(trial,(raw.(n_trial).timepoints{tt}))); 
            
            % Raw Spat. Correlation
            raw.(n_trial).MspatCorr(tt) = ...
                nanmean(MspatCorrTMP(k,raw.(n_trial).timepoints{tt}));
            
            % Raw GEV
            raw.(n_trial).GEV(tt) = ...
                sum( (GFP(trial,raw.(n_trial).timepoints{tt}) .* ...
                MspatCorrTMP(k,raw.(n_trial).timepoints{tt})).^2) ./...
                sum(GFP(trial,:).^2);
            
            % Raw Duration
            raw.(n_trial).Duration(tt) = runs(tt) .* (1000 / fs);
            
        end    
    end
end

%% Transition Probabilities (as with hmmestimate(states,states);)
TP_total = zeros(K);
for trial = 1:Ntrials
    states = order{trial,:};
    states = states(states ~= 0);
    % prepare output matrix
    numStates = K;
    tr = zeros(numStates);
    seqLen = length(states);
    % count up the transitions from the state path
    for count = 1:seqLen-1
        tr(states(count),states(count+1)) = tr(states(count),states(count+1)) + 1;
    end
    trRowSum = sum(tr,2);
    %aggregate transition frequenices across trials
    TP_total = TP_total + tr;
    % if we don't have any values then report zeros instead of NaNs.
    trRowSum(trRowSum == 0) = -inf;
    % normalize to give frequency estimate.
    TP(:,:,trial) = tr./repmat(trRowSum,1,numStates);
end
%divide by rowsum after aggregation
total_RowSum = sum(TP_total, 2);
total_RowSum(total_RowSum == 0) = -inf;
TP_total = TP_total ./ repmat(total_RowSum, 1, K);


%% Write to EEG structure
%   average per trial:
Mstats.GEVtotal = GEVtotal;
Mstats.Gfp = MGFP;
Mstats.Occurence = MOcc;
Mstats.Duration = MDur;
Mstats.Coverage = TCov;
Mstats.GEV = GEV;
Mstats.MspatCorr = MspatCorr;
Mstats.TP = TP;
Mstats.seq = seq;
Mstats.msFirstTF = msFirstTF;
Mstats.polarity = polarity;

% MTruninger added this 17.02.2021
% raw data per trial
Mstats.raw = raw;

if Ntrials > 1
    % mean parameters
    Mstats.avgs.Gfp = nanmean(MGFP,1);
    Mstats.avgs.Occurence = nanmean(MOcc,1);
    Mstats.avgs.Duration = nanmean(MDur,1);
    Mstats.avgs.Coverage = nanmean(TCov,1);
    Mstats.avgs.GEV = nanmean(GEV,1);
    Mstats.avgs.MspatCorr = nanmean(MspatCorr,1);
    Mstats.avgs.TP = TP_total;
    
    % standard deviation of parameters
    Mstats.avgs.stdGfp = nanstd(MGFP);
    Mstats.avgs.stdOccurence = nanstd(MOcc);
    Mstats.avgs.stdDuration = nanstd(MDur);
    Mstats.avgs.stdCoverage = nanstd(TCov);
    Mstats.avgs.stdGEV = nanstd(GEV);
    Mstats.avgs.stdMspatCorr = nanstd(MspatCorr);
end
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

% MTruninger added this function 17.02.2021
function tp = tp_indivMS(r)
%% tp_indivMS
% This function extracts the timepoints of each single microstate based on
% the "runs" (output c from the my_RLE-function).

tp = cell(1,length(r)); %prepare output cell array
curr_t2 = 0; 

for j = 1:length(r) % loop over all single microstates
    t1 = curr_t2 + 1; % get start timepoint of current single microstate
    t2 = t1 + r(j) - 1; % get end timepoint of current single microstate
    curr_t2 = t2; 
    %output
    tp{j} = [t1:t2]; % save all timepoints belonging to the current single microstate
end 
end 
