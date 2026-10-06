%% Appendix1_figure1_source_code1.m
%
% Purpose:
%   Generate individual STED fluorescence intensity line profiles
%   directly from Appendix 1 - Figure 1 source data.
%
% No smoothing or normalization is applied.
%
% =========================================================================

clear;
close all;
clc;


%% ------------------------------------------------------------------------
% USER SETTINGS
% -------------------------------------------------------------------------

input_file = 'Appendix1_Figure1_SourceData1.xlsx';
sheet_name = 1;

% Excel blocks:
% Distance Glass | Glass |
% Distance noFA  | noFA  |
% Distance FA    | FA

profiles(1).name  = 'alpha5 total VC5';
profiles(1).range = 'C:H';
profiles(1).ylim  = [0 8];

profiles(2).name  = '9EG7 active beta1';
profiles(2).range = 'Q:V';
profiles(2).ylim  = [0 18];

profiles(3).name  = 'alphaVbeta3';
profiles(3).range = 'X:AC';
profiles(3).ylim  = [0 18];

profiles(4).name  = 'Paxillin';
profiles(4).range = 'AI:AN';
profiles(4).ylim  = [0 26];


%% ------------------------------------------------------------------------
% COLORS
% -------------------------------------------------------------------------

gray_col = [0.50 0.50 0.50];
noFA_col = [0.85 0.40 0.70];
FA_col   = [0.00 0.45 0.10];


%% ------------------------------------------------------------------------
% LOOP THROUGH THE FOUR PROFILES
% -------------------------------------------------------------------------

for p = 1:numel(profiles)

    fprintf('\n-----------------------------------------\n');
    fprintf('Reading: %s\n', profiles(p).name);
    fprintf('Excel range: %s\n', profiles(p).range);


    %% --------------------------------------------------------------------
    % READ DATA
    % ---------------------------------------------------------------------

    M = readmatrix(input_file, ...
        'Sheet', sheet_name, ...
        'Range', profiles(p).range);


    % Check that six columns were read
    if size(M,2) < 6
        error(['Range %s does not contain 6 columns. ' ...
               'Check the Excel organization.'], ...
               profiles(p).range);
    end


    %% --------------------------------------------------------------------
    % EXTRACT COLUMNS
    % ---------------------------------------------------------------------

    dist_glass = M(:,1);
    glass      = M(:,2);

    dist_noFA  = M(:,3);
    noFA       = M(:,4);

    dist_FA    = M(:,5);
    FA         = M(:,6);


    %% --------------------------------------------------------------------
    % REMOVE EMPTY CELLS / HEADERS
    % ---------------------------------------------------------------------

    idxGlass = isfinite(dist_glass) & isfinite(glass);
    idxNoFA  = isfinite(dist_noFA)  & isfinite(noFA);
    idxFA    = isfinite(dist_FA)    & isfinite(FA);

    dist_glass = dist_glass(idxGlass);
    glass      = glass(idxGlass);

    dist_noFA  = dist_noFA(idxNoFA);
    noFA       = noFA(idxNoFA);

    dist_FA    = dist_FA(idxFA);
    FA         = FA(idxFA);


    %% --------------------------------------------------------------------
    % DIAGNOSTIC
    % ---------------------------------------------------------------------

    fprintf('Glass points = %d\n', numel(glass));
    fprintf('noFA points  = %d\n', numel(noFA));
    fprintf('FA points    = %d\n', numel(FA));


    % Stop if one of the three datasets is empty
    if isempty(glass)
        error('%s: Glass data are empty. Check range %s.', ...
            profiles(p).name, profiles(p).range);
    end

    if isempty(noFA)
        error('%s: noFA data are empty. Check range %s.', ...
            profiles(p).name, profiles(p).range);
    end

    if isempty(FA)
        error('%s: FA data are empty. Check range %s.', ...
            profiles(p).name, profiles(p).range);
    end


    %% --------------------------------------------------------------------
    % CREATE FIGURE AND AXES
    % ---------------------------------------------------------------------

    fig = figure( ...
        'Units','pixels', ...
        'Position',[100 + 50*p, 100 + 30*p, 420, 420], ...
        'Color','w', ...
        'Renderer','painters', ...
        'Name',profiles(p).name);

    ax = axes('Parent',fig);

    hold(ax,'on');


    %% --------------------------------------------------------------------
    % PLOT PROFILES
    % ---------------------------------------------------------------------

    hGlass = plot(ax, ...
        dist_glass, glass, '-', ...
        'Color',gray_col, ...
        'LineWidth',1.7);

    hNoFA = plot(ax, ...
        dist_noFA, noFA, '--', ...
        'Color',noFA_col, ...
        'LineWidth',2.0);

    hFA = plot(ax, ...
        dist_FA, FA, '-', ...
        'Color',FA_col, ...
        'LineWidth',2.4);


    %% --------------------------------------------------------------------
    % AXIS FORMATTING
    % ---------------------------------------------------------------------

    xlim(ax,[0 1]);
    ylim(ax,profiles(p).ylim);

    axis(ax,'square');

    xlabel(ax,'Distance (\mum)', ...
        'FontName','Arial', ...
        'FontSize',18);

    ylabel(ax,'Photon Counts (a.u.)', ...
        'FontName','Arial', ...
        'FontSize',18);

    set(ax, ...
        'FontName','Arial', ...
        'FontSize',14, ...
        'LineWidth',1.2, ...
        'TickDir','out', ...
        'TickLength',[0.02 0.02], ...
        'Box','off');


    %% --------------------------------------------------------------------
    % LEGEND
    % ---------------------------------------------------------------------

    legend(ax, ...
        [hGlass hNoFA hFA], ...
        {'Glass','noFA','FA'}, ...
        'Location','northeast', ...
        'Box','off', ...
        'FontName','Arial', ...
        'FontSize',11);


    %% --------------------------------------------------------------------
    % TITLE
    % ---------------------------------------------------------------------

    title(ax,profiles(p).name, ...
        'FontName','Arial', ...
        'FontSize',12, ...
        'FontWeight','normal');


    %% --------------------------------------------------------------------
    % FIGURE POSITION
    % ---------------------------------------------------------------------

    ax.Units = 'normalized';
    ax.Position = [0.18 0.16 0.76 0.76];


    %% --------------------------------------------------------------------
    % EXPORT
    % ---------------------------------------------------------------------

    export_name = sprintf( ...
        'Appendix1_Figure1C_profile_%d', p);

    exportgraphics(fig, ...
        [export_name '.pdf'], ...
        'ContentType','vector');

    exportgraphics(fig, ...
        [export_name '.png'], ...
        'Resolution',600);

end


fprintf('\n=========================================\n');
fprintf('Finished. Four profiles processed.\n');
fprintf('=========================================\n');