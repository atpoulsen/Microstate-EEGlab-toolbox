% MicroPlotSegments() - plots Microstate segmention over the GFP
%
%  Draws a plot with Microstate segments over the GFP, with a colour for
%  each Microstate. Optionally the microstate numbers over each segment and
%  the topographical microstate maps can by plotted.
%  Note that the GFP plot is not good at handling one microstate
%  appearing for one sample within another (e.g. if a single microstate
%  timepoint is plotted first, it will be "over-plotted" by the surrounding
%  microstate).
%
%  Usage:
%   >> OUTEEG = MicroPlotSegments ( INEEG, 'key1', 'val1', 'key2', 'val2' ... )
%
%  Please cite this toolbox as:
%  Poulsen, A. T., Pedroni, A., Langer, N., &  Hansen, L. K. (unpublished
%  manuscript). Microstate EEGlab toolbox: An introductionary guide.
%
%  Inputs
%  EEG      - EEG-lab EEG structure with the fields 'data' (channels x
%             samples (x epochs)), 'times', 'srate' and 'chanlocs'; and the
%             'data', 'labels' and 'prototypes' fields in EEG.microstates.
%
%  Optional input:
%  'label_type' - Plot labels from 'segmentation' (default) or 'backfit'?.
%  'plotsegnos' - Plot Microstate numbers above microstate segments?
%                 'first' (default) plots a number over the first segment
%                 for each microstate, 'all' plots over all segments.
%                 'none' means no numbers are plotted.
%  'plottopos'  - Plot topography of microstates? 1, yes (default); 0, no.
%  'plot_time'  - [min max] time range, in ms, to be plotted. If empty
%                 (default), all time points will be plotted.
%
% Authors:
% Andreas Pedroni, andreas.pedroni@uzh.ch
% University of Zürich, Psychologisches Institut, Methoden der
% Plastizitätsforschung.
%
% Andreas Trier Poulsen, atpo@dtu.dk
% Technical University of Denmark, DTU Compute, Cognitive systems.
%
% Parts of the code is modified from the timtopo() function by Scott
% Makeig.
%
% April 2017.
%
% See also: timtopo(), topoplot()

% Copyright (C) 2017 Andreas Trier Poulsen, atpo@dtu.dk
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

function MicroPlotSegments(EEG, varargin)
%% Error check and initialisation
if nargin < 1
    help MicroPlotSegments
    return
else
    settings = check_settings(varargin);
end

MAX_TOPOS = 20; % maximum number of topoplots

% check if chanlocs are present
if isempty(EEG.chanlocs)
    error('there are no channel locations, please load them (with pop_readlocs)')
end

% check if segmentation has been run
if ~isfield(EEG.microstate,'data')
    error('No data selected for segmentation. Run "Select data" first.')
end


%% Read data from EEG and EEGlab
chanlocs = EEG.chanlocs;
A = EEG.microstate.prototypes;
K = size(A,2);

% Get microstate labels
if strcmp(settings.label_type, 'segmentation')
    if ~isfield(EEG.microstate, 'labels')
        error(['Label type ''segmentation'' not present in EEG.microstate.' ...
            'Run segmentation or try backfit ''label'' type.'])
    end
    labels = EEG.microstate.labels;
    % load data and times
    if ischar(EEG.microstate.data)
        data = EEG.data;
        times = EEG.times;
    else
        data = EEG.microstate.data;
        times = (1:size(data,2))/EEG.srate*1e3; % Assuming ms.
    end
elseif strcmp(settings.label_type, 'backfit')
    if ~isfield(EEG.microstate.fit, 'labels')
        error(['Backfitted labels not present in EEG.microstate.fit.labels' ...
            'Run backfitting or try labels from segmentation instead.'])
    end
    labels = EEG.microstate.fit.labels;
    data = EEG.data;
    times = EEG.times;
else
    error('Label type ''%s'' not supported.', settings.label_type)
end

if ~isempty(settings.plot_time)
   t_start = find(times>=settings.plot_time(1),1,'first');
   t_end = find(times>=settings.plot_time(2),1,'first');
   
   data = data(:,t_start:t_end);
   times = times(t_start:t_end);
   labels = labels(t_start:t_end);
end

icadefs;

times(end+1) = times(end) + 1e3/EEG.srate; % adding extra sample before the first one. Assuming ms.


%% Compute title and axes font sizes
% figure('Name','Microstates');
pos = get(gca,'Position');
% pos = [0.1300 0.1100 0.7750 0.8150];
axis('off') % don't show axis
cla % clear the current axes
if pos(4)>0.70
    titlefont= 16;
    axfont = 16;
elseif pos(4)>0.40
    titlefont= 14;
    axfont = 14;
elseif pos(4)>0.30
    titlefont= 12;
    axfont = 12;
elseif pos(4)>0.22
    titlefont= 10;
    axfont = 10;
else
    titlefont= 8;
    axfont = 8;
end


%% Compute topoplot head width and separation
if settings.plottopos
    if K > MAX_TOPOS
        fprintf('MicroPlotSegments(): too many Microstates - only first %d will be shown!\n',MAX_TOPOS);
        plottimes = plottimes(1:MAX_TOPOS);
        K = MAX_TOPOS;
    end
    
    head_sep = 0.3;
    topowidth = pos(3)/((6*K-1)/5); % width of each topoplot
    if topowidth> 0.25*pos(4) % dont make too large (more than 1/4 of axes width)!
        topowidth = 0.25*pos(4);
    end
    
    halfn = floor(K/2);
    if rem(K,2) == 1  % odd number of topos
        topoleft = pos(3)/2 - (K/2+halfn*head_sep)*topowidth;
    else % even number of topos
        topoleft = pos(3)/2 - ((halfn)+(halfn-1)*head_sep)*topowidth;
    end
    % topoleft = topoleft - 0.01; % adjust left a bit for colorbar
end


%% Ready GFP plot
if settings.plottopos
    % Place the plot at bottom of the figure and make it's height dependent on topoheight.
    topo_space = .1;
    middle_space = .06;
    pos(4) = 1 - topo_space - topowidth*(1+head_sep) - middle_space - pos(2);
end
axdata = axes('Units', 'Normalized', 'Position', pos, 'FontSize', axfont); % creates an axis
% axdata = axes('Units','Normalized','Position',[pos(1) pos(2) pos(3) 0.6*pos(4)],'FontSize',axfont); % creates an axis

set(axdata,'Color',BACKCOLOR);
limits = get(axdata,'Ylim');
set(axdata,'GridLineStyle',':')
set(axdata,'Xgrid','off')
set(axdata,'Ygrid','on')
axes(axdata)
axcolor = get(gcf,'Color');
set(axdata,'Color',BACKCOLOR);
hold on

% prepare GFP plot
empty_clusters = []; % to keep track of empty microstates
GFP = std(data,[],1);
GFP = [GFP GFP(end)]; % adding extra sample after the last one.

c20 =   [0.368627450980392,0.309803921568627,0.635294117647059;0.232941976907398,0.426964379443341,0.716724269529944;
    0.197103273996820,0.545428702089154,0.740441532379098;0.289164658427885,0.676519464158588,0.684127372366063;
    0.425488941852069,0.775245691220066,0.646374895135081;0.570745493485241,0.831883120082506,0.644760405293702;
    0.712331911655248,0.884123184534485,0.639177467615844;0.852572369568295,0.945354582576408,0.606699084357464;
    0.916245263768755,0.960207509159080,0.609279235239802;0.939109517434401,0.954525564903514,0.691784259410135;
    0.964614344304082,0.940412515225342,0.691500293230307;0.992457455354208,0.899715439687825,0.585852977458640;
    0.995633845898762,0.831741829511967,0.491082885642943;0.993209665846219,0.718091384842371,0.403306245777860;
    0.984967339319457,0.589815920663289,0.324081427293180;0.962493646532822,0.451160196578616,0.265020625541367;
    0.918594329250906,0.348117017024907,0.280748065372818;0.843470509188835,0.253888856898642,0.309426572786864;
    0.746582251164705,0.137135217599700,0.298133552881715;0.619607843137255,0.00392156862745098,0.258823529411765];

if  K > 7 && K < 21
    cmap = c20;
else
    cmap = colormap('lines');
end


%% Plot GFP
fprintf('Plotting GFP segments for each microstate:')
for k = 1:K
    fprintf(' %d',k)
    
    % plot GFP in color
    x = nan(1,length(times)); % adding extra sample before the first one.
    idx = labels(1,:) == k; % finding indices where microstate is active
    if sum(idx)==0
        warning('Microstate cluster %d is has no members. It''s prototype will not be plotted.',k)
        empty_clusters(end+1) = k;
        continue
    end
    % adding one sample in the end of each segment to avoid white holes
    % between segments.
    idx2 = [false idx(1:end-1)] | idx;
    
    % adding extra sample at the end
    if labels(end) == k
        idx2 = [idx2 true];
    else
        idx2 = [idx2 false];
    end
    x(idx2) = GFP(1,idx2);
    
    a = area(times,x);
    a.EdgeColor = cmap(k,:);
    a.FaceColor = cmap(k,:);
    
    % plot microstate numbers
    if sum(strcmp(settings.plotsegnos,{'first','all'}))
        % finding start and end of each segment for MS
        num_gap = max(GFP)*0.025;
        seg_start = find(diff([0 idx])==1);
        seg_end = find(diff([0 idx])==-1);
        if labels(end) == k
            % if last segment belongs to MS
            seg_end(end+1) = length(idx);
        end
        
        % plot number over middle of segment(s)
        if strcmp(settings.plotsegnos,'all')
            for seg = 1:length(seg_start)
                textx = (times(seg_end(seg)) + times(seg_start(seg)))/2;
                seg_ymax = max(GFP((seg_start(seg):seg_end(seg))));
                texty = double(seg_ymax + num_gap);
                %                 texty = double(GFP(num_idx(seg)) + num_gap);
                text(textx, texty, num2str(k), 'FontSize',axfont-3);
            end
        else
            % only first segment
            seg = 1;
            textx = (times(seg_end(seg)) + times(seg_start(seg)))/2;
            seg_ymax = max(GFP((seg_start(seg):seg_end(seg))));
            texty = double(seg_ymax + num_gap);
            %                 texty = double(GFP(num_idx(seg)) + num_gap);
            text(textx, texty, num2str(k), 'FontSize',axfont-3);
        end
        
    end
end
fprintf('.\n')


%% Set labels
xl= xlabel('Latency (ms)');
set(xl,'FontSize',axfont);
xlim([min(times) max(times)])
yl = ylabel('GFP'); %ylabel('GFP (\muV)');
set(yl,'FontSize',axfont,'FontAngle','normal');

% check if xlabel is cut off by bottom of figure
outpos = get(gca,'OuterPosition');
if outpos(2) < 0
    % move axis up, so label fits in figure.
    outpos(2) = 0;
    set(gca,'OuterPosition',outpos);
end


%% Plot the topoplots
if settings.plottopos
    disp('Plotting microstate topoplots')
    topoaxes = zeros(1,K);
    % all microstates, which occur in the segments
    for k = 1:K
        if sum(empty_clusters == k)
            % checking if microstate is empty. If so, skipping topoplot.
            continue
        end
        % [pos(3)*topoleft+pos(1)+(t-1)*(1+head_sep)*topowidth ...
        axtp = axes('Units','Normalized','Position',...
            [topoleft+pos(1)+(k-1)*(1+head_sep)*topowidth, ...
            pos(2)+pos(4)+middle_space, ...
            topowidth, topowidth*(1+head_sep)]);
        axes(axtp) % topoplot axes
        topoaxes(k) = axtp; % save axes handles
        cla
        topoplot(A(:,k),chanlocs, 'electrodes', 'off', 'style', 'map','headrad','rim',...
            'hcolor',cmap(k,:));
        text(0.00,0.80,['MS ', num2str(k)],'FontSize',axfont-3,...
            'HorizontalAlignment','Center'); % ,'fontweight','bold');
    end
end
end

function settings = check_settings(vargs)
%% check settings
% Checks settings given as optional inputs for MicroPlot.
% Undefined inputs is set to default values.
varg_check = {   'label_type'  'string' []  'segmentation' ;
    'plotsegnos'  'string' []  'first' ;
    'plot_time' 'real' [] [] ;
    'plottopos' 'integer' [] 1};
settings = finputcheck( vargs, varg_check);
if ischar(settings), error(settings); end; % check for error
end
