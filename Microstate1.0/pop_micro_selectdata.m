% pop_micro_selectdata() - selects and aggregates data for microstate segmentation.
%
% pop_micro_selectdata selects data for microstate segmentation, in
% severeral different ways:
% * For spontaneous data: selection of GFP peaks.
% * For event related data: ERPs as average over epochs. If aggregating
%   multiple datasets they can either be concatenated or averaged into
%   single ERP. 
% It is possible to aggregate data from multiple datasets loaded in ALLEEG.
% When aggregating data a new dataset 'NewEEG' will be created and added to
% ALLEEG.
%
% Usage:
%   >> EEG = pop_micro_selectdata ( EEG ); % pop up window
%   >> [EEG, ALLEEG, CURRENTSET, NewEEG] = pop_micro_selectdata ( EEG, ...
%         ALLEEG ); % pop up window
%   >> [EEG, ALLEEG, CURRENTSET, NewEEG] = pop_micro_selectdata ( EEG, ...
%         ALLEEG, 'key1', 'val1', 'key2', 'val2' ... )
%
% Please cite this toolbox as:
% Poulsen, A. T., Pedroni, A., Langer, N., &  Hansen, L. K. (unpublished
% manuscript). Microstate EEGlab toolbox: An introductionary guide.
%
% Inputs:
%   EEG              - EEGlab EEG structure.
%
% Optional inputs:
%   ALLEEG           - EEGlab structure containing all datasets read into
%                      EEGlab as EEG structures.
%  'datatype'        - 'spontaneous': Finds GFP peaks. 'ERPavg': Averages 
%                       over epochs and aggregated datasets. 'ERPconc':  
%                       Averages over epochs and concatenates aggregated
%                       datasets. Default is 'ERPconc'.
%   multiple datasets they can either be concatenated or averaged into
%   single ERP.  . -----                      -**********
%  'dataset_idx'     - Indices for datasets in ALLEEG to aggragate data
%                      from. Leave empty to use current dataset (default).
%  'avgref'          - Calculate average reference. 1 - yes (default), 
%                      0 - no.
%  'normalise'       - Normalise each dataset with average channel std. 
%                      1 - yes, 0 - no (default).
%  'Markers'         - String that indicates the marker to which data is
%                      epoched. !Not implemented yet!
%  'MinPeakDist'     - Minimum Distance between GFP peaks in ms
%                      (default = 10 ms).
%  'Npeaks'          - Number of GFP peaks per subject that enter the
%                      segmentation. Note that the maximum number of peaks
%                      is restricted to the minimum number of GFP peaks
%                      across subjects.
%  'GFPthresh'       - Reject peaks over threshold (multiples of std(GFP)).
%                      Set to avoid extreme GFP peaks that are likely 
%                      caused by artifacts. Set to zero to turn off
%                      (default).
% Outputs:
%   EEG        - EEGlab EEG structure with added microstate field.
%   ALLEEG     - EEGlab structure containing all datasets read into
%                EEGlab as EEG structures.
%   CURRENTSET - Workspace variable index of the current dataset. Only
%                relevant when creating a new dataset through aggragating
%                data.
%   NewEEG     - New EEGlab EEG structure containing data selected for
%                microstate analysis and info in added microstate field.
%
% Authors:
% Andreas Trier Poulsen, atpo@dtu.dk
% Technical University of Denmark, DTU Compute, Cognitive systems.
%
% Andreas Pedroni, andreas.pedroni@uzh.ch
% University of Zürich, Psychologisches Institut, Methoden der
% Plastizitätsforschung. 
%
% August 2017.
%
% See also: eeglab

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

function [EEG, ALLEEG, com, NewEEG] = pop_micro_selectdata(EEG, ALLEEG, varargin)
%% Error check and initialisation
com = '';

if nargin < 1
    help pop_micro_selectdata;
    return;
elseif nargin == 1
    ALLEEG = EEG; % in case ALLEEG was not given as input.
end

if nargin < 3
    % pop-up window in case no further input is given
    settings = input_popup(ALLEEG);
    if strcmp(settings,'cancel')
        return
    end
else
    settings = check_settings(varargin);
end


%% Readying settings
dataset_idx = settings.dataset_idx;
if isempty(dataset_idx)
   aggregate_data = 0;
   Ndatasets = 1;
   NewEEG = [];
else
   if length(dataset_idx) == 1
      error('Please select more than one dataset, when aggregating data.') 
   end
   aggregate_data = 1; 
   Ndatasets = length(dataset_idx);
   % check consistency of selected datasets
   NewEEG = check_datasets_consistency(ALLEEG(dataset_idx));
end


%% do the aggregation for spontaneous or erp data
if strcmp(settings.datatype,'spontaneous')
    %% Readying settings for countinuous data 
    Npeaks = settings.Npeaks;
    MinPeakDist  = settings.MinPeakDist;
    GFPthresh = settings.GFPthresh;
    
    % loop through the subjects
    for i=1:Ndatasets
        %% Load data
        if aggregate_data
            X = ALLEEG(dataset_idx(i)).data;
            fs = ALLEEG(dataset_idx(i)).srate;
        else
            X = EEG.data;
            fs = EEG.srate;
        end
        X = reshape(X ,size(X,1),size(X,2)*size(X,3));
        
        
        %% Average reference
        if settings.avgref
            X = bsxfun(@minus, X, mean(X));
        end
        
        
        %% Normalise by average channel std.
        if settings.normalise
            X = X ./ mean(std(X,0,2));
        end
        
        
        %% Parse GFP peaks
        % do the segmentation on the GFP peaks, because they are high in
        % signal-to-noise
        
        % calculate the GFP
        GFP = double(std(X,[],1));
        
        % Minimum distance in ms:
        MinTFdistance = MinPeakDist * fs/1000; % assuming ms
        [~, peakidx{i,1}] = findpeaks(GFP,'MinPeakDistance',MinTFdistance); 
        
        % OPTION, which may be included: to avoid extreme GFP peaks that are likely caused by artifacts
        if GFPthresh > 0
            noisepeak_idx = GFP(peakidx{i}) > (mean(GFP) + GFPthresh*std(GFP));
            peakidx{i} = peakidx{i}(~noisepeak_idx);
        end
    end

    %% Select the GPF peaks 
    % select a number of random peaks (if there are many peaks or you are
    % impatient)
    % this is the maximum number that can be selected of each subject
    maxNpeaks  = min(cellfun('length',peakidx));
    
    if Npeaks > maxNpeaks
        Npeaks = maxNpeaks;
    end
    
    GFPdata = [];
    for i = 1:Ndatasets
        % load data
        if aggregate_data
            X = ALLEEG(dataset_idx(i)).data;
        else
            X = EEG.data;
        end
        X = reshape(X ,size(X,1),size(X,2)*size(X,3));
        
        % Average reference
        if settings.avgref
            X = bsxfun(@minus, X, mean(X));
        end

        % Normalise by average channel std.
        if settings.normalise
            X = X ./ mean(std(X,0,2));
        end
        
        % find a number of random peaks
        selection = randperm(length(peakidx{i}));
        GFPpeakidx{i} = peakidx{i}(selection(1:Npeaks));
        GFPdata = [GFPdata, X(:,GFPpeakidx{i})]; 
    end

    
    %% Save chosen data in the EEG struct
    if aggregate_data
        % Add information to NewEEG struct
        NewEEG.microstate.data = 'set_data';
        NewEEG.microstate.GFPpeakidx = GFPpeakidx;
        NewEEG.microstate.data_origin.ALLEEG_idx = dataset_idx;
        dataset_names = {ALLEEG.setname};
        NewEEG.microstate.data_origin.dataset_names = dataset_names(dataset_idx);
        
        NewEEG.data = GFPdata;
        NewEEG.setname = 'MicroGFPpeakData';
        NewEEG.pnts = size(NewEEG.data,2);
        NewEEG.chanlocs = ALLEEG(dataset_idx(1)).chanlocs; % assuming datasets have the same chanlocs
        NewEEG.times = 0:1000/NewEEG.srate:(length(NewEEG.data)*1000/NewEEG.srate)-1;
        NewEEG.event = [];
        NewEEG.urevent = [];
        NewEEG.eventdescription = [];
        ALLEEG = eeg_store(ALLEEG, NewEEG);
    else
        % Save data in EEG struct given as input.
        EEG.microstate.data = GFPdata;
        EEG.microstate.GFPpeakidx = GFPpeakidx{1};
    end

    
elseif sum(strcmp(settings.datatype,{'ERPavg','ERPconc'}))
    ERP = [];
    for i = 1:Ndatasets
        if aggregate_data
            X = ALLEEG(dataset_idx(i)).data;
        else
            X = EEG.data;
        end

        
        %% Average reference
        if settings.avgref
            X = bsxfun(@minus, X, mean(X));
        end
        
        
        %% Normalise by average channel std.
        if settings.normalise
            X = X ./ mean(std(X(:,:),0,2));
        end
        
        
        set_avg = mean(X,3);
        % preparing datasets in case of aggregation of multiple datasets 
        if strcmp(settings.datatype,'ERPavg')
            % grand average across datasets
            ERP = cat(3, ERP, set_avg);
        else
            % concatenating dataset ERPs
            ERP = cat(2, ERP, set_avg);
        end
    end
    
    if strcmp(settings.datatype,'ERPavg')
        % averaging across datasets in case of aggregation
        ERP = mean(ERP, 3);
    end
    
    
    %% Save chosen data in the EEG struct
    if aggregate_data
        % Add information to NewEEG struct
        NewEEG.microstate.data = 'set_data';
        NewEEG.microstate.data_origin.ALLEEG_idx = dataset_idx;
        dataset_names = {ALLEEG.setname};
        NewEEG.microstate.data_origin.dataset_names = dataset_names(dataset_idx);
        
        NewEEG.data = ERP;
        NewEEG.setname = 'MicroERPdata';
        NewEEG.pnts = size(NewEEG.data,2);
        NewEEG.chanlocs = ALLEEG(dataset_idx(1)).chanlocs; % assuming datasets have the same chanlocs
        NewEEG.times = 0:1000/NewEEG.srate:(length(NewEEG.data)*1000/NewEEG.srate)-1;
        NewEEG.event = [];
        NewEEG.urevent = [];
        NewEEG.eventdescription = [];
        ALLEEG = eeg_store(ALLEEG, NewEEG);
    else
        % Save data in EEG struct given as input.
        EEG.microstate.data = ERP;
    end  
end


%% Define command string
com = sprintf('[%s, %s] = pop_micro_selectdata( %s, %s', inputname(1), ...
    inputname(2),inputname(1), inputname(2));
com = settings_to_string(com,settings);
com = [com ' );'];

end

% -------------- Pop-ups-------------- %
function settings = input_popup(ALLEEG)
% Function for creating popup window to input algorithm settings
%

%% Create Inputs for popup
% Title string
info_str1 = 'Please note that this is an early version of the plugin. Bug-reports and suggestions';
info_str2 = 'are welcome at atpo@dtu.dk.';
line.info = { {'Style' 'text' 'string' info_str1} ...
    {'Style' 'text' 'string' info_str2} {} };
geo.info = {1 1 1};

% Datatype
style.datatype = 'popupmenu';
dropdown_data = {'ERP - Concatenate datasets' 'ERP - Average datasets'  ... % For dropdown menu
    'Spontaneous - GFP peaks'};
popmenu.datatype = {'ERPconc' 'ERPavg' 'spontaneous'}; % Corresponding calls for pop-function
data_str = dropdown_data{1}; %string for popupmenu
for data_idx = 2:length(dropdown_data); data_str = [data_str '|' dropdown_data{data_idx}]; end;
line.datatype = { {'Style' 'text' 'string' 'Data type:'}, ...
    {'Style' style.datatype 'string' data_str 'tag' 'datatype' 'value' 1} };
geo.datatype = {[1 1]};

% Aggregate data?
style.aggregate_data = 'checkbox';
line.aggregate_data = { {'Style' style.aggregate_data 'value' 0 'string' ...
    'Select and aggregate data from other datasets (will open new window).'...
    'tag' 'aggregate_data'} };
geo.aggregate_data = {1};

% Average reference?
style.avgref = 'checkbox';
line.avgref = { {'Style' style.avgref 'value' 1 'string' 'Calculate average reference.' ...
    'tag' 'avgref'} };
geo.avgref = {1};

% Normalise dataset(s)?
style.normalise = 'checkbox';
norm_tipstr = 'Normalise each dataset with average channel std. ';
line.normalise = { {'Style' style.normalise 'value' 0 'string' 'Normalise dataset(s).' ...
    'tooltipstring' norm_tipstr 'tag' 'normalise'} };
geo.normalise = {1};

% GFP info string
gfp_info_str1 = 'Settings for the extraction of GFP peak maps.';
gfp_info_str2 = '(Settings will be ignored if ERP is chosen).';
line.gfp_info = { {} {'Style' 'text' 'string' gfp_info_str1 'fontweight' 'bold'}...
    {'Style' 'text' 'string' gfp_info_str2}};
geo.gfp_info = {1 1 1};

% MinPeakDist
style.MinPeakDist = 'edit';
line.MinPeakDist = { {'Style' 'text' 'string' 'Minimum peak distance (ms):'}, ...
    {'Style' style.MinPeakDist 'string' ' 10 ' 'tag' 'MinPeakDist'} };
geo.MinPeakDist = {[1 .2]};

% Npeaks
style.Npeaks = 'edit';
npeaks_tipstr = ['Note that the maximum number of peaks is restricted to '...
    'the minimum number of GFP peaks across subjects.'];
line.Npeaks = { {'Style' 'text' 'string' ...
    'No. of GFP peaks per subject that enter the segmentation.:' ...
    'tooltipstring' npeaks_tipstr}, ...
    {'Style' style.Npeaks 'string' ' 1000 ' 'tag' 'Npeaks'} };
geo.Npeaks = {[1 .2]};

% GFPthresh
style.GFPthresh = 'edit';
thresh_tipstr = ['Set to avoid extreme GFP peaks that are likely caused by '...
    'artifacts. Set to zero to turn off.'];
line.GFPthresh = { {'Style' 'text' 'string' 'Reject peaks over threshold (multiples of std(GFP)):' ...
    'tooltipstring' thresh_tipstr}, ...
    {'Style' style.GFPthresh 'string' ' 0 ' 'tag' 'GFPthresh'} };
geo.GFPthresh = {[1 .2]};


% Order inputs for GUI
geometry = [geo.info geo.datatype geo.aggregate_data geo.avgref ...
    geo.normalise geo.gfp_info geo.MinPeakDist geo.Npeaks geo.GFPthresh];
uilist = [line.info line.datatype line.aggregate_data line.avgref ...
    line.normalise line.gfp_info line.MinPeakDist line.Npeaks line.GFPthresh];


%% Create Popup
[~,~,~,pop_out] = inputgui( geometry, uilist, ...
    'pophelp(''pop_micro_selectdata'');', ...
    'Select and aggregates data for microstate segmentation -- pop_micro_selectdata()');


%% Interpret output from popup
if isstruct(pop_out)
    settings = struct;
    settings = interpret_popup(pop_out, settings, style, popmenu);
    
    % select the datasets, which should be included in the analysis
    if settings.aggregate_data
        datasetNames = {ALLEEG.setname};
        dataset_idx = listdlg2('ListString',datasetNames,'PromptString',...
            {'Select datasets to aggregate:'});
    else
        dataset_idx = [];
    end
    settings.dataset_idx = dataset_idx;
    settings = rmfield(settings,'aggregate_data');
    
    % remove settings related spontaneous data
    if ~strcmp(settings.datatype,'spontaneous')
       settings = rmfield(settings,'MinPeakDist');
       settings = rmfield(settings,'Npeaks');
       settings = rmfield(settings,'GFPthresh');
    end
else
    settings = 'cancel';
end

end

% -------------- helper functions -------------- %
function settings = interpret_popup(pop_out, settings, style, popmenu)
% Interpret output from pop_up window, "pop_out", and arrange it in
% "settings" struct. The fields in "style" should be the same as in "pop_out"
% (defined as tags in inputgui.m). The struct popmenu is optional and only
% needed if popmenus are used in the pop_up window.

names = fieldnames(style);
for i = 1:length(names)
    switch style.(names{i})
        case 'edit'
            if isempty(pop_out.(names{i})) % empty?
                settings.(names{i}) = [];
            else
                settings.(names{i}) = eval(pop_out.(names{i}));
            end
        case 'checkbox'
            settings.(names{i}) = pop_out.(names{i});
        case 'popupmenu'
            settings.(names{i}) = popmenu.(names{i}){pop_out.(names{i})};
    end
end

end

function settings = check_settings(vargs)
%% check settings
% Checks optional inputs for pop_micro_selectdata() and enters them in
% settings struct. Undefined inputs is set to default values.
varg_check = {  'datatype'         'string'	[]	'ERPconc';
    'avgref'      'integer'	[]	1;
    'normalise' 'integer'	[]	0;
    'dataset_idx'  'real'    []         []};

% Checking if datatype is defined
dat_ind = find(strcmp(vargs, 'datatype'));
if ~isempty(dat_ind)
    if strcmp(lower(vargs{dat_ind+1}), 'spontaneous')
        % Adding check for spontaneous data related settings
        varg_check = [varg_check;
            {'MinPeakDist'  'real'    []         10;
            'Npeaks'  'integer'    []         1000;
            'GFPthresh'  'real'    []         0 } ];
    end
end

settings = finputcheck( vargs, varg_check);
if ischar(settings), error(settings); end; % check for error
end

function NewEEG = check_datasets_consistency(ALLEEG)
%% checks if basic parameters of EEG data are consistent between datasets.
%% and creates new EEG structure that will later be filled depending on the type of data (spontaneous, ERP)

if isequal(ALLEEG(:).nbchan) &&  isequal(ALLEEG(:).srate); % number of channels and sampling rate
    % Create a clean EEG structure that is used to store data.
    NewEEG = eeg_emptyset();
    NewEEG.srate = ALLEEG(1).srate;
    NewEEG.nbchan = ALLEEG(1).nbchan;
else
    error('The datasets differ in the number of channels and/or the sampling rate')  ;
end
end

function com = settings_to_string(com,settings)
% Adds settings struct to existing com string in the form 'key1', 'val1',
% 'key2', 'val2' ... .
% Can handle structs, strings, vectors and scalars. I.e. not matrices.

names = fieldnames(settings);

for i = 1:length(names)
    if isstruct(settings.(names{i})) % struct?
        com = settings_to_string(com,settings.(names{i}));
    elseif isempty(settings.(names{i})) % empty?
        com = [ com sprintf(', ''%s'', []', names{i}) ];
    elseif ischar(settings.(names{i})) % string?
        com = [ com sprintf(', ''%s'', ''%s''', names{i}, settings.(names{i})) ];
    elseif length(settings.(names{i})) > 1 % vector?
        N_elements = length(settings.(names{i}));
        range = max(settings.(names{i})) - min(settings.(names{i})) + 1;
        if  N_elements == range % write vetor as 'min_value:max_value'
            com = [ com sprintf(', ''%s'', %g:%g', names{i}, ...
                min(settings.(names{i})), max(settings.(names{i}))) ];
        else % write vector with individual elements
            com = [ com sprintf(', ''%s'', [%g', names{i}, ...
                settings.(names{i})(1))];
            for n = 2:N_elements
                com = [ com sprintf(',%g', settings.(names{i})(n))];
            end
            com = [ com ']'];
        end
    else % scalar
        com = [ com sprintf(', ''%s'', %g', names{i}, settings.(names{i})) ];
    end
end

end