% pop_micro_segment() - select settings for segmentation into microstates
%
%  Note - Early untested version. (TAAHC not implemented yet).
%
% Function for segmenting EEG into microstates using either Modified
% K-means as descibed in [1], Variational microstates as described in [2],
% ordinary K-means using Matlabs built in function (Stats Toolbox needed)
% or Topographical Atomize-Agglomerate Hierarchical Clustering as described
% in [3,4].
%
% Usage:
%   >> OUTEEG = pop_micro_segment ( INEEG ); % pop up window
%   >> OUTEEG = pop_micro_segment ( INEEG, 'key1', 'val1', 'key2', 'val2' ... )
%
%  Please cite this toolbox as:
%  Poulsen, A. T., Pedroni, A., Langer, N., &  Hansen, L. K. (unpublished
%  manuscript). Microstate EEGlab toolbox: An introductionary guide.
%
% Inputs:
%   INEEG - Input dataset.
%
% Optional inputs:
%  'algorithm'          - String denoting the algorithm used for
%                         segmentation. ['Variational microstate analysis'|
%                         'Modified K-means'|'K-means'|'Topographical
%                         Atomize-Agglomerate Hierarchical Clustering'|]
%                         (default: 'Modified K-means').
%  'Nmicrostates'       - Number of microstates to segment EEG into. If
%                         vector then each value will be tested and the
%                         optimum number of microstates selected based on
%                         measures of fit (default: 3:8).
%  'sorting'            - Method of sorting microstates. ['Global explained
%                         variance'|'Chronological appearance'|'Frequency']
%                         (default: 'Global explained variance').
%  'verbose'            - Print status messages to command window?
%                         (default:1).
%  'normalise'          - Normalise EEG with average channel std.
%                         (default: 1).
%  'temporal_smoothing' - Employ temporal smoothing of microstate labeling?
%                         (default: 0).
%
% Optional algorithm specific inputs:
% * Modified K-means:
%   'smooth_width'   - Temporal smoothing width. Denotes the samples on each side
%                      of current sample (default: 0). Setting to zero
%                      turns smoothing off.
%   'smooth_weight'  - Temporal smoothing weight (default: 5).
%   'Nrepetitions'   - Number of random initialisations of algorithm
%                      (default: 10).
%   'max_iterations' - Maximum number of iterations of algorithm
%                      (default: 1000).
%   'threshold'      - Threshold of convergence based on relative change
%                      in noise variance (default: 1e-6).
%
% * Variational microstates:
%   'sig2_0'         - Prior variance of activations (default: average
%                      EEG channel variance).
%   'p0'             - Probability for having same microstate as last
%                      timepoint (default: 0). Setting to zero turns
%                      smoothing off.
%   'Nrepetitions'   - Number of random initialisations of algorithm
%                      (default: 10).
%   'max_iterations' - Maximum number of iterations of algorithm
%                      (default: 1000).
%   'threshold'      - Threshold of convergence based on relative change
%                      in noise variance (default: 1e-6).
%
% * K-means:
%   'smooth_width'   - Temporal smoothing width. Denotes the samples on each side
%                      of current sample (default: 0). Setting to zero
%                      turns smoothing off.
%   'smooth_weight'  - Temporal smoothing weight (default: 5).
%   'Nrepetitions'   - Number of random initialisations of algorithm
%                      (default: 10).
%   'max_iterations' - Maximum number of iterations of algorithm
%                      (default: 1000).
%   'threshold'      - Threshold of convergence based on relative change
%                      in noise variance (default: 1e-6).
%
% Outputs:
%   OUTEEG - Output dataset. This function saves output in the substruct 
%            OUTEEG.microstate, which contains info general and algorithm-
%            specific settings as well as the results of segemntation. The
%            .Res substruct contains algorithm specific results. See the
%            relevant algorithm function for explanations on variable
%            names.
%
%
%  [1] - Pascual-Marqui, R. D., Michel, C. M., & Lehmann, D. (1995).
%        Segmentation of brain electrical activity into microstates: model
%        estimation and validation. IEEE Transactions on Biomedical
%        Engineering.
%  [2] - (unpublished manuscript). Variational microstate analysis.
%  [3] - Murray, M. M., Brunet, D., & Michel, C. M. (2008). Topographic
%        ERP analyses: A step-by-step tutorial review. Brain Topography.
%  [4] - Brunet, D.(2011). Cartool reference Guide. Cartool v3.51.
%
% Author: Andreas Trier Poulsen, atpo@dtu.dk
% Technical University of Denmark, Cognitive systems - February 2017
%
% See also: eeglab

% If adding new algorithms remember to:
% * Add algorithm specific settings pop_up.
% * Make sure check_settings is updated with new algorithm specific
%   settings.
% * add variable names from Res struct that should be sorted alongside maps
%   and labels to the sort_names_opt cell.
% * Make sure "com" string can handle new settings.

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

function [OUTEEG, com] = pop_micro_segment(EEG, varargin)
%% Error check and initialisation
OUTEEG = EEG;
com = '';

if nargin < 1
    help pop_micro_segment;
    return;
end;


%% pop-up window in case no further input is given
if nargin < 2
    settings = segment_popup();
    if strcmp(settings,'cancel')
        return
    end
else
    settings = check_settings(varargin);
end

%% Datatype
%ERP?


%% Write settings to OUTEEG (overwrites any previous microstate info )
OUTEEG.microstate = settings;


%% Normalise EEG
if settings.normalise
    % normalising by average channel std.
    EEG.data = EEG.data./mean(std(EEG.data,0,2));
end


%% Run selected algorithm
switch settings.algorithm
    case 'Modified K-means'
        % readying algorithm settings
        opts_settings_names = {
            'lambda' 'smooth_weight';
            'b' 'smooth_width';
            'reps' 'Nrepetitions';
            'max_iterations' 'max_iterations';
            'thresh' 'threshold';
            'verbose' 'verbose'};
        opts = algosettings_to_opts(settings,opts_settings_names);
        K_range = settings.algorithm_settings.Nmicrostates;
        
        % running algorithm
        [OUTEEG.microstate.scalp_maps, OUTEEG.microstate.labels, ...
            OUTEEG.microstate.Res] = modkmeans(EEG.data, K_range, opts);
        
        sort_names_opt = {};
        
        
    case 'Variational microstate analysis'
        % readying algorithm settings
        opts_settings_names = {
            'sig2_0' 'sig2_0';
            'p0' 'p0';
            'reps' 'Nrepetitions';
            'max_iterations' 'max_iterations';
            'thresh' 'threshold';
            'verbose' 'verbose'};
        opts = algosettings_to_opts(settings,opts_settings_names);
        K_range = settings.algorithm_settings.Nmicrostates;
        
        % running algorithm
        [OUTEEG.microstate.scalp_maps, OUTEEG.microstate.labels, ...
            OUTEEG.microstate.Res] = varMicro(EEG.data, K_range, opts);
        
        % Res variables that should sorted alongside scalp_maps and labels.
        sort_names_opt = {'S_opt', 'sig2_z_opt'};
        
    case 'K-means'
        % starting K-merans using subfunction
        OUTEEG = run_kmeans(EEG.data,OUTEEG,settings);
        
        % Res variables that should sorted alongside scalp_maps and labels.
        sort_names_opt = {};
        
    otherwise
        error(['selected algorithm,''' settings.algorithm ''' not available'])
end


%% Calculate measures of fit
[KL, W, CV, GEV] = calc_fitmeas(EEG.data, OUTEEG.microstate.Res.A_all, ...
    OUTEEG.microstate.Res.L_all);

OUTEEG.microstate.Res.KL = KL;
OUTEEG.microstate.Res.W = W;
OUTEEG.microstate.Res.CV = CV;
OUTEEG.microstate.Res.GEV = GEV;


%% Sorting microstates according to chosen method
OUTEEG = sort_microstates(EEG.data,OUTEEG,sort_names_opt);


%% Define command string
com = sprintf('%s = pop_micro_segment( %s', inputname(1), inputname(1));
com = settings_to_string(com,settings);
com = [com ' );'];

% if settings.algorithm_settings.verbose
%    disp('Done.') 
% end

end

% ------------------------------ Pop-ups ---------------------------------%
function settings = segment_popup()
% Function for creating popup window to input algorithm settings
%

%% Create Inputs for popup
% Title string
info_str1 = 'Please note that this is an early version of the plugin. Bug-reports and suggestions';
info_str2 = 'are welcome at atpo@dtu.dk.';
info_str3 = 'This early version assumes continous data i.e, that EEG.data is 2D.';
line.info = { {'Style' 'text' 'string' info_str1} ...
    {'Style' 'text' 'string' info_str2}...
    {'Style' 'text' 'string' info_str3} {} };
geo.info = {1 1 1 1};


% Choose algorithm
style.algorithm = 'popupmenu';
popmenu.algorithm = {'K-means' 'Modified K-means'  ...
    'Topographical Atomize-Agglomerate Hierarchical Clustering' ...
    'Variational microstate analysis'};
algo_str = popmenu.algorithm{1}; %string for popupmenu
for a = 2:length(popmenu.algorithm); algo_str = [algo_str '|' popmenu.algorithm{a}]; end;
line.algorithm = { {'Style' 'text' 'string' 'Choose algorithm:'}, ...
    {'Style' style.algorithm 'string' algo_str 'tag' 'algorithm' 'value' 2} };
geo.algorithm = {[1 1]};

% K range (Nmicrostates)
style.Nmicrostates = 'edit';
line.Nmicrostates = { {'Style' 'text' 'string' 'Number of microstates:'}, ...
    {'Style' style.Nmicrostates 'string' ' 3:8 ' 'tag' 'Nmicrostates'} {},... %end of first line
    {'Style' 'text' 'string' '(If vector, algorithm is run for each value).'},...
    {} }; %end of second line
geo.Nmicrostates = {[1 .3 .7] [1 1]};

% Sorting
style.sorting = 'popupmenu';
popmenu.sorting = {'Global explained variance' 'Chronological appearance' ...
    'Frequency'};
sort_str = popmenu.sorting{1}; %string for popupmenu
for s = 2:length(popmenu.sorting); sort_str = [sort_str '|' popmenu.sorting{s}]; end;
line.sorting = { {'Style' 'text' 'string' 'Sort microstates according to ... :'}, ...
    {'Style' style.sorting 'string' sort_str 'tag' 'sorting'} };
geo.sorting = {[1 1]};

% Smoothing
style.temporal_smoothing = 'checkbox';
smooth_tipstr = ['Employ temporal smoothing of microstate labeling. If '...
    'selected a popup window will ask for smoothing parameters.'];
line.temporal_smoothing = { {'Style' style.temporal_smoothing 'value' 0 'string' 'Temporal smoothing' ...
    'tooltipstring' smooth_tipstr 'tag' 'temporal_smoothing'} {} };
geo.temporal_smoothing = {[1 1]};

% Normalisation
style.normalise = 'checkbox';
norm_tipstr = ['Normalise EEG with average channel std. (recommended '...
    'especially for variational microstates).'];
line.normalise = { {'Style' style.normalise 'value' 1 'string' 'Normalise EEG before segmentation' ...
    'tooltipstring' norm_tipstr 'tag' 'normalise'} {} };
geo.normalise = {[1 1]};

% Verbose
style.verbose = 'checkbox';
line.verbose = { {'Style' style.verbose 'value' 1 'string' ...
    'Print status in command window?' 'tag' 'verbose'} {} };
geo.verbose = {[1 1]};


%% Order inputs for GUI
geometry = [geo.info geo.algorithm geo.sorting {1} geo.Nmicrostates ...
    geo.temporal_smoothing geo.normalise geo.verbose];
uilist = [line.info line.algorithm line.sorting {{}} line.Nmicrostates ...
    line.temporal_smoothing line.normalise line.verbose];


%% Create Popup
[~,~,~,pop_out] = inputgui( geometry, uilist, ...
    'pophelp(''pop_micro_segment'');', 'Segement into microstates -- pop_micro_segment()');


%% Interpret output from popup
if isstruct(pop_out)
    settings = struct;
    settings = interpret_popup(pop_out, settings, style, popmenu);
    % Moving relevant fields to algorithm_settings substruct
    to_algoset = {'Nmicrostates', 'verbose'};
    for i = 1:length(to_algoset)
        settings.algorithm_settings.(to_algoset{i}) = settings.(to_algoset{i});
        settings = rmfield(settings,(to_algoset{i}));
    end
else
    settings = 'cancel';
end


%% Create popup for algorithm specific settings
if ~strcmp(settings,'cancel')
    switch settings.algorithm
        case 'Modified K-means'
            settings = modk_popup(settings);
        case 'Variational microstate analysis'
            settings = varmicro_popup(settings);
        case 'K-means'
            settings = kmeans_popup(settings);
        otherwise
            error(['selected algorithm,''' settings.algorithm ''' not available'])
    end
end

end

function settings = modk_popup(settings)
% Popup for unique input for the modified K-means algorithm
%

%% Create Inputs for popup
% Title string
info_str = 'Input parameters specific for ''Modified K-means''.';
line.info = { {'Style' 'text' 'string' info_str} {} };
geo.info = {1 1};

% Repetitions
style.Nrepetitions = 'edit';
line.Nrepetitions = { {'Style' 'text' 'string' 'No. of random initialisations:'}, ...
    {'Style' style.Nrepetitions 'string' ' 10 ' 'tag' 'Nrepetitions'} };
geo.Nrepetitions = {[1 .2]};

% Max iterations
style.max_iterations = 'edit';
line.max_iterations = { {'Style' 'text' 'string' 'Max. no. of iterations:'}, ...
    {'Style' style.max_iterations 'string' ' 1000 ' 'tag' 'max_iterations'} };
geo.max_iterations = {[1 .2]};

% Threshold
style.threshold = 'edit';
line.threshold = { {'Style' 'text' 'string' 'Relative threshold for convergence:'}, ...
    {'Style' style.threshold 'string' ' 1e-6 ' 'tag' 'threshold'} };
geo.threshold = {[1 .2]};

% Smoothing info string
smooth_info_str = 'Temporal smoothing.';
line.smooth_info = { {} {'Style' 'text' 'string' smooth_info_str 'fontweight' 'bold'}};
geo.smooth_info = {1 1};

% Smoothing width
style.smooth_width = 'edit';
if_width = fastif(settings.temporal_smoothing, { ' 3 ' }, { '0' 'enable' 'off' });
line.smooth_width = { {'Style' 'text' 'string' 'Width of smoothing windows (in samples):'}, ...
    {'Style' style.smooth_width 'string' if_width{:} 'tag' 'smooth_width'} };
geo.smooth_width = {[1 .2]};

% Smoothing weight
style.smooth_weight = 'edit';
if_weight = fastif(settings.temporal_smoothing, { ' 5 ' }, { '0' 'enable' 'off' });
line.smooth_weight = { {'Style' 'text' 'string' 'Smoothing weight:'}, ...
    {'Style' style.smooth_weight 'string' if_weight{:} 'tag' 'smooth_weight'} };
geo.smooth_weight = {[1 .2]};


%% Order inputs for GUI
geometry = [geo.info geo.Nrepetitions geo.max_iterations geo.threshold ...
    geo.smooth_info geo.smooth_width geo.smooth_weight];
uilist = [line.info line.Nrepetitions line.max_iterations line.threshold ...
    line.smooth_info line.smooth_width line.smooth_weight];


%% Create Popup
[~,~,~,pop_out] = inputgui( geometry, uilist, ...
    'pophelp(''pop_micro_segment'');', 'Extra input for modified K-means -- pop_micro_segment()');


%% Interpret output from popup
if isstruct(pop_out)
    settings.algorithm_settings = interpret_popup(pop_out, settings.algorithm_settings, style);
else
    settings = 'cancel';
end
end

function settings = varmicro_popup(settings)
% Popup for unique input for the variatinal microstates algorithm
%
%

%% Create Inputs for popup
g = [1 .3];
% Title string
info_str = 'Input parameters specific for ''Variational Microstates''.';
line.info = { {'Style' 'text' 'string' info_str} {} };
geo.info = {1 1};

% Repetitions
style.Nrepetitions = 'edit';
line.reps = { {'Style' 'text' 'string' 'No. of random initialisations:'}, ...
    {'Style' style.Nrepetitions 'string' ' 10 ' 'tag' 'Nrepetitions'} };
geo.reps = {g};

% Max iterations
style.max_iterations = 'edit';
line.max_iter = { {'Style' 'text' 'string' 'Max. no. of iterations:'}, ...
    {'Style' style.max_iterations 'string' ' 1000 ' 'tag' 'max_iterations'} };
geo.max_iter = {g};

% Threshold
style.threshold = 'edit';
line.thresh = { {'Style' 'text' 'string' 'Relative threshold for convergence:'}, ...
    {'Style' style.threshold 'string' ' 1e-6 ' 'tag' 'threshold'} };
geo.thresh = {g};

% Prior variance of activations (sig2_z0)
style.sig2_0 = 'edit';
line.sig2_0 = { {'Style' 'text' 'string' 'Prior variance of activations (sig2_z0): '}, ...
    {'Style' style.sig2_0 'string' '' 'tag' 'sig2_0'},... %end of first line
    {'Style' 'text' 'string' '(Positive value or leave empty for data variance).'},...
    {} }; %end of second line
geo.sig2_0 = {g g};

% Smoothing info string
smooth_info_str = 'Temporal smoothing';
line.smooth_info = { {} {'Style' 'text' 'string' smooth_info_str 'fontweight' 'bold'}};
geo.smooth_info = {1 1};

% Smoothing weight
style.p0 = 'edit';
if_weight = fastif(settings.temporal_smoothing, { ' 0.3 ' }, { '0' 'enable' 'off' });
line.weight = { {'Style' 'text' 'string' 'Probability for having same microstate as last timepoint (p0):'}, ...
    {'Style' style.p0 'string' if_weight{:} 'tag' 'p0'},... %end of first line
    {'Style' 'text' 'string' '(smoothing weight, between 0 and 1).'},...
    {}  };
geo.weight = {g g};


%% Order inputs for GUI
geometry = [geo.info geo.reps geo.max_iter geo.thresh geo.sig2_0 geo.smooth_info geo.weight];
uilist = [line.info line.reps line.max_iter line.thresh line.sig2_0 line.smooth_info line.weight];


%% Create Popup
[~,~,~,pop_out] = inputgui( geometry, uilist, ...
    'pophelp(''pop_micro_segment'');', 'Extra input for variational microstates -- pop_micro_segment()');


%% Interpret output from popup
if isstruct(pop_out)
    settings.algorithm_settings = interpret_popup(pop_out, settings.algorithm_settings, style);
else
    settings = 'cancel';
end


end

function settings = kmeans_popup(settings)
% Popup for unique input for the modified K-means algorithm
%

%% Create Inputs for popup
% Title string
info_str1 = 'Input parameters specific for ''K-means''.';
info_str2 = 'Note that classic ''K-means'' is not polarity-invariant. ''K-means'' is therefore likely';
info_str3 = 'to require twice as many scalp maps to explain the EEG.';
line.info = { {'Style' 'text' 'string' info_str1} ...
    {'Style' 'text' 'string' info_str2} ...
    {'Style' 'text' 'string' info_str3} {} };
geo.info = {1 1 1 1};

% Repetitions
style.Nrepetitions = 'edit';
line.Nrepetitions = { {'Style' 'text' 'string' 'No. of random initialisations:'}, ...
    {'Style' style.Nrepetitions 'string' ' 10 ' 'tag' 'Nrepetitions'} };
geo.Nrepetitions = {[1 .2]};

% Max iterations
style.max_iterations = 'edit';
line.max_iterations = { {'Style' 'text' 'string' 'Max. no. of iterations:'}, ...
    {'Style' style.max_iterations 'string' ' 1000 ' 'tag' 'max_iterations'} };
geo.max_iterations = {[1 .2]};

% Smoothing info string
smooth_info_str1 = 'Temporal smoothing.';
smooth_info_str2 = 'Note: For K-means the mean squared error is minimised directly.';
line.smooth_info = { {} {'Style' 'text' 'string' smooth_info_str1 'fontweight' 'bold'}, ...
    {'Style' 'text' 'string' smooth_info_str2} };
geo.smooth_info = {1 1 1};

% Smoothing width
style.smooth_width = 'edit';
if_width = fastif(settings.temporal_smoothing, { ' 3 ' }, { '0' 'enable' 'off' });
line.smooth_width = { {'Style' 'text' 'string' 'Width of smoothing windows (in samples):'}, ...
    {'Style' style.smooth_width 'string' if_width{:} 'tag' 'smooth_width'} };
geo.smooth_width = {[1 .2]};

% Smoothing weight
style.smooth_weight = 'edit';
if_weight = fastif(settings.temporal_smoothing, { ' 5 ' }, { '0' 'enable' 'off' });
line.smooth_weight = { {'Style' 'text' 'string' 'Smoothing weight:'}, ...
    {'Style' style.smooth_weight 'string' if_weight{:} 'tag' 'smooth_weight'} };
geo.smooth_weight = {[1 .2]};

% Threshold
style.threshold = 'edit';
if_thresh = fastif(settings.temporal_smoothing, { ' 1e-6 ' }, { 'NaN' 'enable' 'off' });
line.threshold = { {'Style' 'text' 'string' 'Relative threshold for convergence:'}, ...
    {'Style' style.threshold 'string' if_thresh{:} 'tag' 'threshold'} };
geo.threshold = {[1 .2]};


%% Order inputs for GUI
geometry = [geo.info geo.Nrepetitions geo.max_iterations ...
    geo.smooth_info geo.smooth_width geo.smooth_weight geo.threshold];
uilist = [line.info line.Nrepetitions line.max_iterations ...
    line.smooth_info line.smooth_width line.smooth_weight line.threshold];


%% Create Popup
[~,~,~,pop_out] = inputgui( geometry, uilist, ...
    'pophelp(''pop_micro_segment'');', 'Extra input for K-means -- pop_micro_segment()');


%% Interpret output from popup
if isstruct(pop_out)
    settings.algorithm_settings = interpret_popup(pop_out, settings.algorithm_settings, style);
else
    settings = 'cancel';
end
end
% ----------------------------------------------------------------------- %

% -------------------------- Helper functions --------------------------- %
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
% Checks settings given as optional inputs for pop_micro_segment.
% The function checks and rearranges optional inputs to pop_micro_segment.m
% struct. Undefined inputs is set to default values.
%% General inputs
varg_check = { 'algorithm'  'string'    []         'Modified K-means';
    'Nmicrostates'  'real'    []         3:8;
    'sorting'  'string'    []         'Global explained variance';
    'verbose'  'integer'    []         1;
    'normalise' 'integer'    []         1
    'temporal_smoothing' 'integer'    []         0};

% Fields to be moved to algorithm_settings substruct
to_algoset = {'Nmicrostates', 'verbose'};

%% Algorithm specifik inputs
algo = find(strcmp('algorithm',vargs)) + 1;
switch vargs{algo}
    case 'Modified K-means'
        varg_check = [varg_check;
            { 'smooth_weight'  'real'    []         [];
            'smooth_width'  'real'    []         [];
            'Nrepetitions'  'integer'    []         10;
            'max_iterations'  'integer'    []         1000;
            'threshold'  'real'    []         1e-6 } ];
        to_algoset = [ to_algoset {'smooth_weight', 'smooth_width', ...
            'Nrepetitions', 'max_iterations', 'threshold'} ];
    case 'Variational microstate analysis'
        varg_check = [varg_check;
            { 'Nrepetitions'  'integer'    []         10;
            'max_iterations'  'integer'    []         1000;
            'threshold'  'real'    []         1e-6;
            'sig2_0'  'real'    []         [];
            'p0'  'real'    []         [] } ];
        to_algoset = [ to_algoset {'Nrepetitions', 'max_iterations', ...
            'threshold', 'sig2_0', 'p0'} ];
    case 'K-means'
        varg_check = [varg_check;
            { 'smooth_weight'  'real'    []         [];
            'smooth_width'  'real'    []         [];
            'Nrepetitions'  'integer'    []         10;
            'max_iterations'  'integer'    []         1000;
            'threshold'  'real'    []         1e-6 } ];
        to_algoset = [ to_algoset {'smooth_weight', 'smooth_width', ...
            'Nrepetitions', 'max_iterations', 'threshold'} ];
    otherwise
        error(['selected algorithm,''' vargs{algo} ''' not available'])
end
settings = finputcheck( vargs, varg_check);
if ischar(settings), error(settings); end; % check for error


%% Moving relevant fields to algorithm_settings substruct
for i = 1:length(to_algoset)
    settings.algorithm_settings.(to_algoset{i}) = settings.(to_algoset{i});
    settings = rmfield(settings,(to_algoset{i}));
end

end

function opts = algosettings_to_opts(settings,names)
% Arranges and renames relevant algorithm settings into opts struct.
% "names" is a Nx2 cell, where each row is a opts.name1 = settings.name2
% pair. Does not carry empty fields over

for i = 1:size(names,1)
    if ~isempty(settings.algorithm_settings.(names{i,2})) % checking for empty fields
        opts.(names{i,1}) = settings.algorithm_settings.(names{i,2});
    end
end
end

function OUTEEG = run_kmeans(X,OUTEEG,settings)
% Runs K-means clustering using Matlabs standard implementation.
%

%% Renaming settings and getting ready
K_range = settings.algorithm_settings.Nmicrostates;
reps = settings.algorithm_settings.Nrepetitions;
verbose = settings.algorithm_settings.verbose;
max_iter = settings.algorithm_settings.max_iterations;

MSE_mcv_opt = inf;
K_ind = 0;
if verbose
    disp('Starting K-means')
end

%% Preallocating
N_K = length(K_range);
A_all = cell(N_K,1);
L_all = cell(N_K,1);
MSE = nan(N_K,1);
MSE_mcv = nan(N_K,1);


%% Looping over Number of microstates
for K = K_range
    K_ind = K_ind + 1;
    if verbose
        fprintf(['Analysis no. %i out of %i. Starting %i random '...
            'initialisations for %i microstates.\n'], ...
            K_ind,length(K_range),reps, K)
    end
    
    % Segmenting with ordinary K-means clustering
    if verbose
        [L, A] = kmeans(X', K, 'MaxIter', max_iter, 'Replicates', reps,...
            'Display','final');
    else
        [L, A] = kmeans(X', K, 'MaxIter', max_iter, 'Replicates', reps);
    end
    L=L'; A=A';
    
    if settings.temporal_smoothing
        L = smoothing(X,A,L,K,settings);
    end
    
    
    A_all{K_ind} = A;
    L_all{K_ind} = L;
    
    % Calculating 
    [MSE(K_ind), MSE_mcv(K_ind)] = calc_MSE(X,A,L,K);
    
    % Saving optimum solution amongst different values of K
    if MSE_mcv(K_ind) < MSE_mcv_opt
        OUTEEG.microstate.scalp_maps = A;
        OUTEEG.microstate.labels = L;
        MSE_mcv_opt = MSE_mcv(K_ind);
        K_act = K;
    end
    
end


%% Saving to Res struct
OUTEEG.microstate.Res.K_act = K_act;
OUTEEG.microstate.Res.A_all = A_all;
OUTEEG.microstate.Res.L_all = L_all;
OUTEEG.microstate.Res.MSE = MSE;
OUTEEG.microstate.Res.MSE_mcv = MSE_mcv;


end

function [MSE, MSE_mcv] = calc_MSE(X,A,L,K)
% Calculating measures of fit for Kmeans. 

% const = sum(sum(X.^2));
[C,N] = size(X);

% Calculating modified version of the predictive residual variance
% presented in [1] eq. (22).
MSE = mean(mean((X-A(:,L)).^2));
MSE_mcv = MSE * ((C-1)^-1 * (C-1-K))^-2;

% % Noise variance as calculated in Pascual-Marqui (1995)
% sig2_modk = (const - sum(sum(A(:,L).*X).^2)) / (N*(C-1));
% sig2_modk_mcv = sig2_modk * ((C-1)^-1 * (C-1-K))^-2;
% sig2_D = const / (N*(C-1));
% R2 = 1 - sig2_modk/sig2_D;


end

function L = smoothing(X,A,L,K,settings)
%  Modified version of the Segmentation Smoothing Algorithm, as described
%  in Table II of [1]. This version seeks to directly minimise the mean
%  squared error between each EEG sample and the microstate assigned to
%  that timepoint. Smoothing is done using the interval t-b to t+b
%  excluding t.
%  Note, that temporary allocation of labels (denoted with Lambda in [1])
%  is not necessary in this implementation, and steps 3 and 6 are therefore
%  left out.

%% Initialisation (step 1 to 4)
% Reading settings
lambda = settings.algorithm_settings.smooth_weight;
b = settings.algorithm_settings.smooth_width;
max_iterations = settings.algorithm_settings.max_iterations;
thresh = settings.algorithm_settings.threshold;
verbose = settings.algorithm_settings.verbose;

if verbose
    fprintf('Smoothing... ')
end

[C,N] = size(X);

% Check to avoid the loop getting caught by switching one label back and
% forth between iterations.
L_old{1} = zeros(size(L));
L_old{2} = zeros(size(L));

% Step 1
MSE_mod_old = 0;
MSE_mod = Inf;

% Step 4
e = sum(sum((X-A(:,L)).^2)) / (N*(C-1));

% Defining constant for step 5b (replacing original term with MSE)
const_5b = nan(K,N);
for k = 1:K
    const_5b(k,:) = sum(bsxfun(@minus,X,A(:,k)).^2) / (2*e*(C-1));
end


%% Iterations (step 5 to 8)
ind = 0;
while abs(MSE_mod_old-MSE_mod) >= thresh*MSE_mod && max_iterations>ind ...
        && mean(L_old{rem(ind,2)+1} == L)~=1
    ind = ind + 1;
    MSE_mod_old = MSE_mod;
    L_old{abs(rem(ind,2)-2)} = L;
    
    % Step 5a
    Nbkt_tmp = zeros(K,N);
    for k = 1:K
        Nbkt_tmp(k,:) = double(L==k);
    end
    % Using filter to count the number of labels equal to k before (tmp1)
    % and after (tmp2) a given timepoint.
    tmp1 = filter([0 ones(1,b)],1,Nbkt_tmp,[],2);
    tmp2 = filter([0 ones(1,b)],1,Nbkt_tmp(:,end:-1:1),[],2);
    Nbkt = tmp1 + tmp2(:,end:-1:1);
    
    % Step 5b
    [~,L] = min( const_5b - lambda*Nbkt );
    
    % Step 7
    MSE_mod = sum(sum((X-A(:,L)).^2)) / (N*(C-1));
end

%% Done
if verbose
    fprintf('Smoothing done in %i iterations.\n',ind)
end

end

function OUTEEG = sort_microstates(X, OUTEEG, sort_names_opt)
% Sorting microstates according to chosen method. GFP and GEV is
% implemented as descibed in [3]. Sorting Z_all, A_all and L_all (from Res.)
% as default, looping over all K in K_range. scalp_maps, labels are also
% sorted by assigning the relevant variables from A_ll{K_act} and
% L_all{K_act}. Variables defined in sort_names_opt are sorted for K_act 
% (variables needs to have K as its first dimension).

if OUTEEG.microstate.algorithm_settings.verbose
    fprintf('Sorting microstates using %s \n',OUTEEG.microstate.sorting)
end

%% Reading settings
K_range = OUTEEG.microstate.algorithm_settings.Nmicrostates;


%% Looping over K_range
K_ind = 0;
for K = K_range
    K_ind = K_ind + 1;
    
    % Assigning A and L
    A = OUTEEG.microstate.Res.A_all{K_ind};
    L = OUTEEG.microstate.Res.L_all{K_ind};
    
    % check to see if a solution was found, otherwise skipping
    if isempty(A)
        continue
    end
    
    %% Calculating chosen sorting measure
    switch OUTEEG.microstate.sorting
        case 'Global explained variance'
            meas_name = 'GEVk';
            sort_method = 'descend';
            
            GFP = var(X);
            % Calculating and summing GEV over all timepoints for the
            % corresponding active microstates
            C2 = columncorr(X,A(:,L));
            GEV = (GFP.*C2).^2;
            sortmeas = zeros(K,1);
            for k=1:K
                sortmeas(k) = sum(GEV(L==k))/sum(GFP.^2);
            end
            
        case 'Chronological appearance'
            meas_name = 'First_appearance';
            sort_method = 'ascend';
            [~,sortmeas] = unique(L);
            
        case 'Frequency'
            meas_name = 'Microstate_frequency';
            sort_method = 'descend';
            sortmeas = nan(K,1);
            for k = 1:K
                sortmeas(k) = sum(L==k);
            end
            
        otherwise
            error('something went wrong during sorting of microstates.')
    end
    
    
    %% Rearranging according to sorting measure
    [sortmeas,idx] = sort(sortmeas,sort_method);
    
    % Saving to sorting measure to OUTEEG
    OUTEEG.microstate.Res.(meas_name){K_ind} = sortmeas;
    
    % Scalp maps
    OUTEEG.microstate.Res.A_all{K_ind} = A(:,idx);
    % Labels
    for k = 1:K
        OUTEEG.microstate.Res.L_all{K_ind}(L==idx(k)) = k;
    end
    
    % Z (not for ordinary Kmeans)
    if ~strcmp(OUTEEG.microstate.algorithm,'K-means')
        OUTEEG.microstate.Res.Z_all{K_ind} = OUTEEG.microstate.Res.Z_all{K_ind}(idx,:);
    end
    
    
    %% Sorting for variables only available for K_act.
    if K == OUTEEG.microstate.Res.K_act
        % scalp_maps and labels
        OUTEEG.microstate.scalp_maps = OUTEEG.microstate.Res.A_all{K_ind};
        OUTEEG.microstate.labels = OUTEEG.microstate.Res.L_all{K_ind};
        
        % Optional sorting for selected Res variables. NOTE! the Res.(variable)
        % needs to have K as its first dimension.
        for i = 1:length(sort_names_opt)
            OUTEEG.microstate.Res.(sort_names_opt{i}) = ...
                OUTEEG.microstate.Res.(sort_names_opt{i})(idx,:);
        end
    end
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