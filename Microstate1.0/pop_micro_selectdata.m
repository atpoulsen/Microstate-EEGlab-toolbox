% pop_micro_selectdata() - selects and aggregates data for microstate segmentation.
%
%  Note - Early untested version.
%
% MicroAggr aggregates data for microstate segmentation. For spontaneous
% data: selection of GFP peaks. For Event related data: grand averaging. 
%
% Usage:
%   >> EEG = pop_micro_selectdata ( EEG ); % pop up window
%   >> [EEG, ALLEEG] = pop_micro_selectdata ( EEG, ALLEEG ); % pop up window
%   >> [EEG, ALLEEG] = pop_micro_selectdata ( EEG, ALLEEG, 'key1', 'val1', 'key2', 'val2' ... )
%
% Please cite this toolbox as:
% Poulsen, A. T., Pedroni, A., Langer, N., &  Hansen, L. K. (unpublished
% manuscript). Microstate EEGlab toolbox: An introductionary guide.
%
% Inputs:
%   EEG             - EEGlab EEG structure.
%
% Optional inputs:
%   ALLEEG           - EEGlab structure containing all datasets read into
%                      EEGlab as EEG structures.
%  'datatype'        - 'Continuous','ERP'.
%  'dataset_idx'     - Indices for datasets in ALLEEG to aggragate data
%                      from. Leave empty to use current dataset (default).
%  'avgref'          - Calculate average reference. 1 - yes (default), 
%                      0 - no.
%  'normalise'       - Normalise each dataset with average channel std. 
%                      1 - yes (default), 0 - no.
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
%   EEG    - xxxxxxxxxxxx
%   ALLEEG - xxxxxxxxx
%
% Authors:
% Andreas Trier Poulsen, atpo@dtu.dk
% Technical University of Denmark, DTU Compute, Cognitive systems.
%
% Andreas Pedroni, andreas.pedroni@uzh.ch
% University of Zürich, Psychologisches Institut, Methoden der
% Plastizitätsforschung. 
%
% May 2017.
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

function [EEG, ALLEEG, com] = pop_micro_selectdata(EEG, ALLEEG, varargin)
%% Error check and initialisation
com = '';

if nargin < 1
    help pop_micro_selectdata;
    return;
elseif nargin == 1
    ALLEEG = EEG; % in case ALLEEG was not given as input.
end;

if nargin < 3
    % pop-up window in case no further input is given
    settings = input_popup(ALLEEG);
    if strcmp(settings,'cancel')
        return
    end
else
    settings = check_settings(varargin);
end;


%% Readying settings
dataset_idx = settings.dataset_idx;
Npeaks = settings.Npeaks;
MinPeakDist  = settings.MinPeakDist;
GFPthresh = settings.GFPthresh;

if isempty(dataset_idx)
   aggregate_data = 0;
else
   aggregate_data = 1; 
end


%% do the aggregation for continuous or erp data
if strcmp(settings.datatype,'Continuous')
    
    % loop through the subjects
    for i=1:length(dataset_idx)
        %% Load data
        if aggregate_data
            X = ALLEEG(dataset_idx(i)).data;
        else
            X = EEG.data;
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
%         plot(GFP)
        
        %% OPTION, which may be included: to avoid extreme GFP peaks that are likely caused by artifacts
        if GFPthresh > 0
            GFP = GFP(GFP < GFPthresh*std(GFP));
        end
        % Minimum distance in tf:
        MinTFdistance = MinPeakDist * ALLEEG(dataset_idx(i)).srate/1000; 
        [~, peakidx{i,1}] = findpeaks(GFP,'MinPeakDistance',MinTFdistance); 
        % check if the peak searching algorithm does something appropriate
        % findpeaks(GFP(1:500),'MinPeakDistance',1,'Annotate','extents');
    end

    %% Select the GPF peaks 
    % select a number of random peaks (if there are many peaks or you are
    % impatient)
    % this is the maximum number that can be selected of each subject
    maxNpeaks  = min(cellfun('length',peakidx));
    
    if ~exist('numSelPeaks','var')
        Npeaks = maxNpeaks;
    else
        if Npeaks > maxNpeaks
            Npeaks = maxNpeaks;
        end
    end
    
    GFPdata = [];
    for i=1:length(dataset_idx)
        % load data
        if aggregate_data
            X = ALLEEG(dataset_idx(i)).data;
        else
            X = EEG.data;
        end
        X = reshape(X ,size(X,1),size(X,2)*size(X,3));
        
        
        % find a number of random peaks
        selection = randperm(length(peakidx{i,1}));
        GFPpeaks(i,:) = peakidx{1,1}(1,selection(1:Npeaks));
        GFPdata = [GFPdata, X(GFPpeaks(i,:))]; 
    end

    
    %% Save chosen data in the EEG struct
    if aggregate_data
        % Create a clean EEG structure that is used to store data.
        NewEEG = struct;
        NewEEG.microstate.data = 'EEGdata';
        NewEEG.microstate.GFPpeaks = GFPpeaks(:)';
        NewEEG.microstate.data_origin.ALLEEG_idx = dataset_idx;
        dataset_names = {ALLEEG.setname};
        NewEEG.microstate.data_origin.dataset_names = dataset_names(dataset_idx);
        
        NewEEG.data = GFPdata;
        NewEEG.setname = 'MicroGFPpeakData';
        NewEEG.pnts = size(NewEEG.data,2);
        NewEEG.nbchans = 1;
        NewEEG.times = 0:1000/NewEEG.srate:(length(NewEEG.data)*1000/NewEEG.srate)-1;
        NewEEG.event = [];
        NewEEG.urevent = [];
        NewEEG.eventdescription = [];
        [ALLEEG, EEG] = eeg_store(ALLEEG, NewEEG);
    else
        % Save data in EEG struct given as input.
        EEG.microstate.data = GFPdata;
        EEG.microstate.GFPpeaks = GFPpeaks(:)';
    end

    
elseif strcmp(settings.datatype,'ERP')
    GA = [];
    for i=1:length(dataset_idx)
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
            X = X ./ mean(std(X,0,2));
        end
        
        % calculate the grand average
        GA = cat(3,GA,X);
    end
    
    
    %% Save chosen data in the EEG struct
    if aggregate_data
        % Create a clean EEG structure that is used to store data.
        NewEEG = struct;
        NewEEG.microstate.data = 'EEGdata';
        NewEEG.microstate.data_origin.ALLEEG_idx = dataset_idx;
        dataset_names = {ALLEEG.setname};
        NewEEG.microstate.data_origin.dataset_names = dataset_names(dataset_idx);
        
        NewEEG.data = mean(GA,3);
        NewEEG.setname = 'MicroERPdata';
        NewEEG.pnts = size(NewEEG.data,2);
        NewEEG.nbchans = 1;
        NewEEG.times = 0:1000/NewEEG.srate:(length(NewEEG.data)*1000/NewEEG.srate)-1;
        NewEEG.event = [];
        NewEEG.urevent = [];
        NewEEG.eventdescription = [];
        [ALLEEG, EEG] = eeg_store(ALLEEG, NewEEG);
    else
        % Save data in EEG struct given as input.
        EEG.microstate.data = mean(GA,3);
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
popmenu.datatype = {'ERP' 'Continuous'};
data_str = popmenu.datatype{1}; %string for popupmenu
for data_idx = 2:length(popmenu.datatype); data_str = [data_str '|' popmenu.datatype{data_idx}]; end;
line.datatype = { {'Style' 'text' 'string' 'Datatype:'}, ...
    {'Style' style.datatype 'string' data_str 'tag' 'datatype' 'value' 1} };
geo.datatype = {[1 .3]};

% Aggregate data?
style.aggregate_data = 'checkbox';
line.aggregate_data = { {'Style' style.aggregate_data 'value' 0 'string' ...
    'Select and aggregated data from other datasets (will open new window).'...
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
line.normalise = { {'Style' style.normalise 'value' 1 'string' 'Normalise dataset(s).' ...
    'tooltipstring' norm_tipstr 'tag' 'normalise'} };
geo.normalise = {1};

% GFP info string
gfp_info_str = 'Settings for continuous data (will be ignored if ERP is chosen).';
line.gfp_info = { {} {'Style' 'text' 'string' gfp_info_str 'fontweight' 'bold'}};
geo.gfp_info = {1 1};

% MinPeakDist
style.MinPeakDist = 'edit';
line.MinPeakDist = { {'Style' 'text' 'string' 'Minimum peak distance (ms):'}, ...
    {'Style' style.MinPeakDist 'string' ' 10 ' 'tag' 'MinPeakDist'} };
geo.MinPeakDist = {[1 .2]};

% Npeaks
style.Npeaks = 'edit';
npeaks_tipstr = ['Note that the maximum number of peaks is restricted to '...
    'the minimum number of GFP peaks across subjects.'];
line.Npeaks = { {'Style' 'text' 'string' 'Relative threshold for convergence:'}, ...
    {'Style' style.Npeaks 'string' ' 1000 ' 'tooltipstring' npeaks_tipstr ...
    'tag' 'Npeaks'} };
geo.Npeaks = {[1 .2]};

% GFPthresh
style.GFPthresh = 'edit';
thresh_tipstr = ['Set to avoid extreme GFP peaks that are likely caused by '...
    'artifacts. Set to zero to turn off.'];
line.GFPthresh = { {'Style' 'text' 'string' 'Reject peaks over threshold (multiples of std(GFP)):'}, ...
    {'Style' style.GFPthresh 'string' ' 0 ' 'tooltipstring' thresh_tipstr ... 
    'tag' 'GFPthresh'} };
geo.GFPthresh = {[1 .2]};


%% Order inputs for GUI
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
            'Select datasets to aggregate:');
    else
        dataset_idx = [];
    end
    settings.dataset_idx = dataset_idx;
    settings = rmfield(settings,'aggregate_data');
    
    % remove settings related continuous data
    if ~strcmp(settings.datatype,'Continuous')
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
% Checks settings given as optional inputs for MicroFit.
% Undefined inputs is set to default values.
varg_check = {  'datatype'         'char'	[]	'ERP';
    'avgref'      'integer'	[]	1;
    'normalise' 'integer'	[]	1;
    'dataset_idx'  'real'    []         []};

if ~isfield(settings, 'datatype'), settings.datatype = 'ERP'; end
if strcmp(settings.datatype,'continuous');
    % Adding check for continuous data related settings
    varg_check = [varg_check;
        {'MinPeakDist'  'real'    []         10;
        'Npeaks'  'integer'    []         1000;
        'GFPthresh'  'real'    []         0 } ];
end

settings = finputcheck( vargs, varg_check);
if ischar(settings), error(settings); end; % check for error
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