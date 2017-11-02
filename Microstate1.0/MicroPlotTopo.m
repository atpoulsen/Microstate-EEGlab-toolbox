% MicroPlotTopo() - plots microstate topographies for a range of microstate numbers.
%
%  Draws a plot with topographical maps for a range of microstate
%  prototypes numbers.
%
%  Usage:
%   >> OUTEEG = MicroPlotSegments ( INEEG, 'key1', 'val1', 'key2', 'val2' ... )
%
%  Please cite this toolbox as:
%  Poulsen, A. T., Pedroni, A., Langer, N., &  Hansen, L. K. (unpublished
%  manuscript). Microstate EEGlab toolbox: An introductionary guide.
%
%  Inputs
%  EEG      - EEG-lab EEG structure (channels x samples (x epochs)) with
%             the field 'chanlocs'; and 'prototypes' fields in
%             EEG.microstates.
%
%  Optional input:
%  'plot_range'  - The range of numbers of microstates to be displayed (as
%                  vector). If empty (default), measures for all microstate
%                  numbers will be plotted.
%
% Author:
% Andreas Trier Poulsen, atpo@dtu.dk
% Technical University of Denmark, DTU Compute, Cognitive systems.
%
% Parts of the code is modified from the timtopo() function by Scott
% Makeig.
%
% August 2017.
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

function h_topo = MicroPlotTopo(EEG, varargin)
%% Error check and initialisation
if nargin < 1
    help MicroPlotTopos
    return
else
    settings = check_settings(varargin);
end;

MAX_TOPO_WIDTH = 16; % maximum number of microstates to be plotted

% check if chanlocs are present
if isempty(EEG.chanlocs)
    error('there are no channel locations, please load them (with pop_readlocs)')
end

% Check plot_range
if isempty(settings.plot_range)
    plot_range = EEG.microstate.algorithm_settings.Nmicrostates;
else
    plot_range = settings.plot_range;
end

[~, plot_idx] = ismember(plot_range, EEG.microstate.algorithm_settings.Nmicrostates);
if sum(plot_idx==0)
    error('Not all microstate numbers in plot_range is present in Nmicrostates.')
end

if max(plot_range) > MAX_TOPO_WIDTH
    warning('MicroPlotTopos(): too many Microstates to plot in a row - will only plot up to %d microstates.\n', MAX_TOPO_WIDTH);
    plot_range = plot_range( plot_range <= MAX_TOPO_WIDTH );
end


%% Read data from EEG and EEGlab
chanlocs = EEG.chanlocs;
A_all = EEG.microstate.Res.A_all(plot_idx);
K_max = max(plot_range);
Nrows = length(plot_range);
icadefs;


%% Compute title and axes font sizes
% figure('Name','Microstates');
pos = get(gca,'Position');
% pos = [0.1300 0.1100 0.7750 0.8150];
axis('off') % don't show axis
cla % clear the current axes
if pos(4)>0.70
    axfont = 16;
    ylab = 'Number of microstates';
elseif pos(4)>0.40
    axfont = 14;
    ylab = 'Number of microstates';
elseif pos(4)>0.30
    axfont = 12;
    ylab = 'Number of microstates';
elseif pos(4)>0.22
    axfont = 10;
    ylab = 'No. of micros';
else
    axfont = 8;
    ylab = 'No. of micros';
end


%% Compute topoplot head width and separation
head_sep_x = 0.1; % x-axis distance between topos as a ratio of topo_width
head_sep_y = 0.3; % y-axis distance between topos as a ratio of topo_width
ylabel_space = pos(3)/(K_max*3); % space for y axis labels
% calculate width/height of topoplots
topowidth = (pos(3) - ylabel_space) / K_max * (1-head_sep_x*(K_max-1)/K_max); % width of each topoplot
topoheight = pos(4)/Nrows * (1-head_sep_y*(Nrows-1)/Nrows); % height of each topoplot

% topowidth = min(topowidth,topoheight); % set to the smallest size.

% recompute head seperations in the x and y direcetions
x_sep = (pos(3) - ylabel_space - K_max*topowidth) / (K_max-1);

if Nrows~=1
    y_sep = (pos(4) - Nrows*topoheight) / (Nrows-1);
else
    y_sep = 0;
end


%% Plot the topoplots
disp('Plotting microstate topoplots')
h_topo = nan(Nrows,K_max);

% Plot ylabel using text
axylab = axes('Units','Normalized','Position',...
    [pos(1), pos(2), ylabel_space/2, pos(4)]);
axes(axylab)
axis off;
text(0, .5, ylab,'FontSize',axfont,...
    'Rotation',90,'fontweight','bold','HorizontalAlignment','Center');

% Loop over rows and collumns
for r = 1:Nrows
    a_ind = Nrows-r+1; % Starting with last in plot_range
    A = A_all{a_ind}; 
    
    % plot microstate number in y-axis
    % microstate label
    axy = axes('Units','Normalized','Position',...
        [pos(1)+ylabel_space/2, ...
        pos(2)+pos(4) - (topoheight*r + (r-1)*y_sep), ...
        ylabel_space/2, topoheight]);
    axes(axy)
    axis off;
    text(0, .5, num2str(plot_range(a_ind)),'FontSize',axfont,...
        'VerticalAlignment','Middle','fontweight','bold');
    
    % plot topos for microstate number
    for k = 1:plot_range(a_ind)
        %     if sum(empty_clusters == k)
        %         % checking if microstate is empty. If so, marking in text colour.
        %         topo_colour = 'r';
        %         warning('Some of the clusters are microstates don')
        %     else
        %         topo_colour = 'k';
        %     end
        
        axtp = axes('Units','Normalized','Position',...
            [ylabel_space+pos(1) + (k-1)*(x_sep+topowidth), ...
            pos(2)+pos(4) - (topoheight*r + (r-1)*y_sep), ...
            topowidth, topoheight]);
        axes(axtp) % topoplot axes
        h_topo(r,k) = axtp; % save axes handles
        cla
        
        % topoplot
        topoplot(A(:,k),chanlocs, 'electrodes', 'off', 'style', 'map');
        
        % microstate label
        text(0.00,0.5,num2str(k),'FontSize',axfont*0.6,...
            'HorizontalAlignment','Center','VerticalAlignment','Bottom',...
            'fontweight','bold');
        
        drawnow
    end
end
end

function settings = check_settings(vargs)
%% check settings
% Checks settings given as optional inputs for MicroPlot.
% Undefined inputs is set to default values.
varg_check = {   'plot_range' 'real' [] []};
settings = finputcheck( vargs, varg_check);
if ischar(settings), error(settings); end; % check for error
end
