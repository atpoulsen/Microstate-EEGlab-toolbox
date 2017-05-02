% [OUTEEG] = MicroAggr(EEG,Datatype,varargin) 

% MicroAggr aggregates data for microstate segmentation
% for Spontaneous data: selection of GFP peaks
% for Event related data: grand averaging 

% 'Datatype' - 'Continuous','ERP'
% 'Markers' - String that indicates the marker to which data is epoched
% 'MinPeakDistance' - Minimum Distance between GFP peaks in ms (default = 10 ms)
% 'numSelPeaks' - Number of GFP peaks per subject that enter the
%                 segmentation. Note that the maximum number of peaks is restricted to the
%                 minimum number of GFP peaks across subjects. 


function [EEG, ALLEEG] = MicroAggr(EEG, Datatype, ALLEEG, varargin)

%% Error check and initialisation
if nargin < 3
    help MicroAggr;
    return;
end;


% select the datasets, which should be included in the analysis
if length(ALLEEG) > 1
datasetNames = {ALLEEG.setname};
Selection = listdlg2('ListString',datasetNames,'PromptString','Select files of one condition:');
else 
Selection = 1;
end

% defaults:
MinPeakDistance = 10;
numSelPeaks = 1000
avgref = 1;
normalise = 1;


%% do the aggregation for continuous or erp data
if strcmp(Datatype,'Continuous')
    % loop through the subjects
    for i=1:length(Selection)
        %% Load data
        X = ALLEEG(Selection(i)).data;
        X  = reshape(X ,size(X,1),size(X,2)*size(X,3));
        
        
        %% Average reference
        if avgref
            X = bsxfun(@minus, X, mean(X));
        end
        
        
        %% Normalise by average channel std.
        if normalise
            EEG.data = EEG.data./mean(std(EEG.data,0,2));
        end
        
        
        %% Parse GFP peaks
        % do the segmentation on the GFP peaks, because they are high in
        % signal-to-noise
        
        % calculate the GFP
        GFP = double(std(X,[],1));
        plot(GFP)
        
        %% OPTION, which may be included: to avoid extreme GFP peaks that are likely caused by artifacts
        % GFP = GFP(GFP < 4*std(GFP));
        
        % Minimum distance in tf:
        MinTFdistance = MinPeakDistance * ALLEEG(Selection(i)).srate/1000; 
        [~, peakidx{i,1}] = findpeaks(GFP,'MinPeakDistance',MinTFdistance); 
        % check if the peak searching algorithm does something appropriate
        % findpeaks(GFP(1:500),'MinPeakDistance',1,'Annotate','extents');
    end

    %% Select the GPF peaks 
    % select a number of random peaks (if there are many peaks or you are
    % impatient)
    % this is the maximum number that can be selected of each subject
    maxSelPeaks  = min(cellfun('length',peakidx));
    
    if ~exist('numSelPeaks','var')
        numSelPeaks = maxSelPeaks;
    else
        if numSelPeaks > maxSelPeaks
            numSelPeaks = maxSelPeaks;
        end
    end
    
    GFPdata = [];
    for i=1:length(Selection)
        X = ALLEEG(Selection(i)).data;
        selection = randperm(length(peakidx{i,1}));
        GFPpeaks(i,:) = peakidx{1,1}(1,selection(1:numSelPeaks));
        warning('X is from the last dataset only')
        GFPdata = [GFPdata, X(GFPpeaks(i,:))]; 
    end
    
    % create a clean EEG structure that is used as input for the
    % Microstate segmentation
    EEG.data = GFPdata;
    EEG.microstates.GFPpeaks = GFPpeaks(:)';
    EEG.setname = 'MicroGFPpeakData';
    EEG.pnts = size(EEG.data,2);
    EEG.nbchans = 1;
    EEG.times = 0:1000/EEG.srate:(length(EEG.data)*1000/EEG.srate)-1;
    EEG.event = [];
    EEG.urevent = [];
    EEG.eventdescription = [];
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    

elseif strcmp(Datatype,'ERP')
    GA = [];
    for i=1:length(Selection)
        
        X = ALLEEG(Selection(i)).data;
        %% Average reference
        if avgref
            X = bsxfun(@minus, X, mean(X));
        end
        
        
        %% Normalise by average channel std.
        if normalise
            EEG.data = EEG.data./mean(std(EEG.data,0,2));
        end
        % calculate the grand average 
        GA = cat(3,GA,X);
    end
    
    
    %% create a clean EEG structure that is used as input for the Microstate segmentation
    EEG.data = mean(GA,3);
    EEG.setname = 'MicroERPdata';
    EEG.event = [];
    EEG.urevent = [];
    EEG.eventdescription = [];
    [ALLEEG OUTEEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    
end

end

% -------------- helper functions -------------- %

function settings = check_settings(vargs)
%% check settings
% Checks settings given as optional inputs for MicroFit.
% Undefined inputs is set to default values.
varg_check = {  'DataType'         'char'	[]	0;
    'polarity'      'integer'	[]	0;
    'sequentialize' 'integer'	[]	0;
    'getorder'      'integer'	[]	1};

settings = finputcheck( vargs, varg_check);
if ischar(settings), error(settings); end; % check for error
end

