%% =========================================================================
% DNA-PAINT DBSCAN cluster analysis
%
% Purpose:
%   Identification and quantification of DNA-PAINT localization clusters
%   within focal adhesion (FA) and non-FA regions.
%
% Input:
%   Binary localization files (.bin) generated from DNA-PAINT analysis.
%
% Analysis:
%   - localization coordinates are extracted from each ROI;
%   - DBSCAN clustering is applied to XY localization coordinates;
%   - cluster-level and ROI-level metrics are calculated.
%
% DBSCAN parameters:
%   Pixel size: 160 nm
%   Epsilon:    0.1
%   MinPts:     10
%
% Output:
%   - Cluster-level statistics
%   - ROI-level summary CSV files
%   - Global summary tables
%

% =========================================================================

clear all; close all; clc;

% ---------------- Paths & archivos ----------------
addpath(genpath('\\nas01\\SMB02\\Amaris\\3-Analysis tools\\Codes for amaris\\Matlab_lib'))

% >>> Ajusta tu carpeta de ROI <<<
pathname    = 'V:\Amaris\6-DNAPaint_Data_Analyzed\24h on FN\12G10\FA';
filePattern = fullfile(pathname, '*_allChs.bin');   % patrón de archivos ROI
theFiles    = dir(filePattern);

% ---------------- Parámetros ----------------
pixel_size_nm     = 160;     % nm/px
% DBSCAN (en px): usar ~1–1.5x precisión de localización (en px)
EPSILON_px        = 0.10;    % 0.15 px ~ 24 nm con 160 nm/px
MinPts            = 10;      % ajusta según densidad (8–15 típico)

% ROI (multi-región por rasterización robusta)
roi_area_method      = 'mask-multi';   % 'mask-multi' | 'convex'
mask_scale           = 3;              % 2–5 típico (sobremuestreo de malla)
mask_dilate_px       = 0;            % 0–2 px (suaviza bordes)
mask_close_px        = 1;            % 0–3 px (cierra huecos pequeños)
min_region_area_um2  = 0.002;          % filtra islitas

% Bins para modas (histos)
binw_locs                 = 5;
binw_diameter_nm          = 10;
binw_area_um2             = 0.01;
binw_density_locs_per_um2 = 50;

% Resumen global
summaryTable = {};

% =======================================================================
for f = 1:length(theFiles)
    filename = theFiles(f).name;
    fpath    = fullfile(pathname, filename);
    fprintf('Procesando ROI (sin colapso): %s\n', filename);

    % -------- Lectura de localizaciones --------
    X = [];
    try
        params = FindClustersStruct();
        params.use_channels = 1;
        params.i3file = fpath;
        X = get_locs_from_bin(params);   % espera al menos [x y], a veces [x y frame ...]
    catch
        warning('get_locs_from_bin falló en %s. Intentando lectura [x y] float32...', filename);
        try
            X = read_float_pairs_xy(fpath);  % solo [x y]
        catch
            warning('No se pudo parsear %s. Saltando.', filename);
            continue;
        end
    end

    % Validación mínima
    if isempty(X) || size(X,2) < 2
        warning('Archivo %s sin columnas suficientes. Saltando.', filename);
        continue;
    end

    % -------- SIN COLAPSO: usar XY tal cual --------
    collapsed = X(:,1:2);   % NO colapsar

    % -------- DBSCAN en XY --------
    plotif = 0;
    ClustData = DetermineClustersDBSCAN(collapsed, EPSILON_px, MinPts, plotif);
    % ClustData: [ClusterID, x_px, y_px]  (ClusterID==0 => noise)

    % -------- ROI por rasterización (clusters + noise) --------
    all_pts = collapsed;
    [ROI_regions, ROI_area_nm2, ROI_area_um2, ROI_sampler] = ...
        roi_mask_multi(all_pts, pixel_size_nm, roi_area_method, mask_scale, ...
                       mask_dilate_px, mask_close_px, min_region_area_um2);

    % -------- Visualización --------
    fig1 = figure('Visible','on'); hold on;
    if ~isempty(ROI_regions)
        for ir = 1:numel(ROI_regions)
            plot(ROI_regions(ir).X, ROI_regions(ir).Y, 'r-', 'LineWidth', 1.0);
        end
    end
    ClusterIDs = unique(ClustData(:,1));
    colors = lines(max(1,length(ClusterIDs)));
    for iID = 1:length(ClusterIDs)
        cid = ClusterIDs(iID);
        pts = ClustData(ClustData(:,1) == cid, 2:3);
        if cid == 0
            plot(pts(:,1), pts(:,2), '.', 'Color', [0.6 0.6 0.6], 'MarkerSize', 3); % noise
        else
            plot(pts(:,1), pts(:,2), '.', 'Color', colors(iID,:), 'MarkerSize', 6);
        end
    end
    axis image; set(gca, 'YDir', 'reverse');
    xlabel('X (px)'); ylabel('Y (px)');
    title(['Clusters (ROI sin colapso): ', filename], 'Interpreter', 'none');

    % -------- Métricas por clúster --------
    Centroids = [];  % [ClusterID, Cx_px, Cy_px, Area_nm2, NumLocs]
    if any(ClustData(:,1)>0)
        for iID = 1:length(ClusterIDs)
            cid = ClusterIDs(iID);
            if cid==0, continue; end
            pts = ClustData(ClustData(:,1) == cid, 2:3);
            if size(pts,1) >= 3
                try
                    k  = boundary(pts(:,1), pts(:,2));
                    xv = pts(k,1); yv = pts(k,2);
                    rep = [false; (diff(xv)==0 & diff(yv)==0)];
                    xv(rep) = []; yv(rep) = [];
                    if numel(xv) >= 3
                        A_nm2 = polyarea(xv, yv) * (pixel_size_nm^2);
                        cx = mean(xv); cy = mean(yv);
                        plot(cx, cy, '+', 'Color', colors(iID,:), 'MarkerSize', 10, 'LineWidth', 1.2);
                        Centroids = [Centroids; cid, cx, cy, A_nm2, size(pts,1)];
                    end
                catch
                    % ignora hulls problemáticos
                end
            end
        end
    end

    % -------- Métricas a nivel archivo --------
    totalLocs = size(all_pts,1);
    noiseLocs = sum(ClustData(:,1)==0);
    N_clusters = numel(unique(ClustData(ClustData(:,1)>0, 1)));
    ClusterDensity_per_um2 = N_clusters / max(ROI_area_um2, eps);

    % Texto en figura
    txt = sprintf('ROI area: %.3f \\mum^2 | N_{clusters}: %d | Cluster density: %.2f cl/\\mum^2', ...
                   ROI_area_um2, N_clusters, ClusterDensity_per_um2);
    annotation('textbox',[0.12 0.01 0.8 0.06],'String',txt, ...
               'EdgeColor','none','Interpreter','tex','FontSize',9);

    % -------- Estadística por clúster + histos --------
    if ~isempty(Centroids)
        LocsPerCluster = Centroids(:,5);
        Areas_nm2      = Centroids(:,4);
        Areas_um2      = Areas_nm2 / 1e6;
        Diameters_nm   = 2 * sqrt(Areas_nm2 / pi);
        Radii_nm       = Diameters_nm / 2;
        Density_lu     = LocsPerCluster .* 1e6 ./ max(Areas_nm2, eps);  % locs/um^2

        meanLocs   = mean(LocsPerCluster);
        medianLocs = median(LocsPerCluster);
        modeLocs   = mode(LocsPerCluster);

        [modeDia, ~]       = modal_bin_center(Diameters_nm,   binw_diameter_nm);
        meanDia            = mean(Diameters_nm);
        medianDia          = median(Diameters_nm);

        [modeArea_um2, ~]  = modal_bin_center(Areas_um2,      binw_area_um2);
        meanArea_um2       = mean(Areas_um2);
        medianArea_um2     = median(Areas_um2);

        [modeDen, ~]       = modal_bin_center(Density_lu,     binw_density_locs_per_um2);
        meanDen            = mean(Density_lu);
        medianDen          = median(Density_lu);

        fig2 = figure('Visible','off','Position',[100,100,1800,700]);
        subplot(2,2,1)
        histogram(LocsPerCluster, 'BinWidth', binw_locs, 'Normalization', 'probability')
        title({['# Localizaciones por clúster'], ...
               ['Mean=', num2str(meanLocs,'%.1f'), ', Median=', num2str(medianLocs,'%.1f'), ', Mode=', num2str(modeLocs)]})
        xlabel('# locs'); ylabel('Frecuencia relativa');

        subplot(2,2,2)
        histogram(Diameters_nm, 'BinWidth', binw_diameter_nm, 'Normalization', 'probability')
        title({['Diámetro equivalente (nm)'], ...
               ['Mean=', num2str(meanDia,'%.1f'), ', Median=', num2str(medianDia,'%.1f'), ', Mode=', num2str(modeDia,'%.1f')]} )
        xlabel('Diámetro [nm]'); ylabel('Frecuencia relativa');

        subplot(2,2,3)
        histogram(Areas_um2, 'BinWidth', binw_area_um2, 'Normalization', 'probability')
        title({['Área (μm^2)'], ...
               ['Mean=', num2str(meanArea_um2,'%.3f'), ', Median=', num2str(medianArea_um2,'%.3f'), ', Mode=', num2str(modeArea_um2,'%.3f')]} )
        xlabel('Área [μm^2]'); ylabel('Frecuencia relativa');

        subplot(2,2,4)
        histogram(Density_lu, 'BinWidth', binw_density_locs_per_um2, 'Normalization', 'probability')
        title({['Densidad (locs/μm^2)'], ...
               ['Mean=', num2str(meanDen,'%.1f'), ', Median=', num2str(medianDen,'%.1f'), ', Mode=', num2str(modeDen,'%.1f')]} )
        xlabel('locs/μm^2'); ylabel('Frecuencia relativa');

        statsPath = fullfile(pathname, [filename(1:end-4), '_Stats.png']);
        exportgraphics(fig2, statsPath, 'Resolution', 300);
        close(fig2);

        EquivDiameter_nm      = Diameters_nm;
        EquivRadius_nm        = Radii_nm;
        Density_locs_per_um2  = Density_lu;
    else
        LocsPerCluster = []; Areas_nm2 = []; Areas_um2 = [];
        Diameters_nm = []; EquivDiameter_nm = []; EquivRadius_nm = [];
        Density_locs_per_um2 = [];
        meanLocs=NaN; medianLocs=NaN; modeLocs=NaN;
        meanDia=NaN; medianDia=NaN; modeDia=NaN;
        meanArea_um2=NaN; medianArea_um2=NaN; modeArea_um2=NaN;
        meanDen=NaN; medianDen=NaN; modeDen=NaN;
        warning('Sin clústeres válidos en: %s', filename);
    end

    % -------- Guardado --------
    clusterFigPath = fullfile(pathname, [filename(1:end-4), '_Clusters.fig']);
    savefig(fig1, clusterFigPath);
    close(fig1);

    save(fullfile(pathname, [filename(1:end-4), '_ClusterStats.mat']), ...
         'Centroids','LocsPerCluster','Diameters_nm','Areas_nm2', ...
         'EquivDiameter_nm','EquivRadius_nm','Density_locs_per_um2', ...
         'ROI_area_nm2','ROI_area_um2','N_clusters', ...
         'roi_area_method','mask_scale','mask_dilate_px','mask_close_px','min_region_area_um2', ...
         'ClusterDensity_per_um2','ROI_regions','ROI_sampler','-v7.3');

    if ~isempty(Centroids)
        T = array2table([Centroids, EquivDiameter_nm, EquivRadius_nm, Density_locs_per_um2, Areas_um2], 'VariableNames', ...
            {'ClusterID','CentroidX_px','CentroidY_px','Area_nm2','NumLocs', ...
             'EquivDiameter_nm','EquivRadius_nm','Density_locs_per_um2','Area_um2'});
    else
        T = array2table(zeros(0,9), 'VariableNames', ...
            {'ClusterID','CentroidX_px','CentroidY_px','Area_nm2','NumLocs', ...
             'EquivDiameter_nm','EquivRadius_nm','Density_locs_per_um2','Area_um2'});
    end
    writetable(T, fullfile(pathname, [filename(1:end-4), '_ClusterStats.csv']));

    FileMetrics = table(ROI_area_um2, N_clusters, ClusterDensity_per_um2, ...
                        totalLocs, noiseLocs, ...
                        'VariableNames', {'ROI_Area_um2','N_clusters','ClusterDensity_per_um2', ...
                                          'TotalLocs','NoiseLocs'});
    writetable(FileMetrics, fullfile(pathname, [filename(1:end-4), '_FileMetrics.csv']));

    summaryTable = [summaryTable; {
        filename, totalLocs, noiseLocs, ...
        meanLocs,   medianLocs,   modeLocs, ...
        meanDia,    medianDia,    modeDia, ...
        meanArea_um2, medianArea_um2, modeArea_um2, ...
        meanDen,    medianDen,    modeDen, ...
        ROI_area_um2, N_clusters, ClusterDensity_per_um2
    }];
end

% -------- Resumen global --------
summaryPath = fullfile(pathname, 'Resumen_Global.csv');
if exist(summaryPath, 'file'), delete(summaryPath); end
summaryTable = cell2table(summaryTable, 'VariableNames', ...
    {'Filename','TotalLocs','NoiseLocs', ...
     'MeanLocsPerCluster','MedianLocsPerCluster','ModeLocsPerCluster', ...
     'MeanDiameter_nm','MedianDiameter_nm','ModeDiameter_nm', ...
     'MeanArea_um2','MedianArea_um2','ModeArea_um2', ...
     'MeanDensity_locs_per_um2','MedianDensity_locs_per_um2','ModeDensity_locs_per_um2', ...
     'ROI_Area_um2','NumClusters','ClusterDensity_per_um2'});
writetable(summaryTable, summaryPath);
fprintf('Guardado resumen global: %s\n', summaryPath);

% ======================= Helpers =======================

function [ROI_regions, ROI_area_nm2, ROI_area_um2, ROI_sampler] = ...
         roi_mask_multi(all_pts, pixel_size_nm, roi_area_method, mask_scale, ...
                        mask_dilate_px, mask_close_px, min_region_area_um2)

    ROI_regions = struct('X',{},'Y',{});
    ROI_area_nm2 = NaN; ROI_area_um2 = NaN;
    ROI_sampler = struct('type','','bbox',[NaN NaN NaN NaN],'scale',NaN);

    if isempty(all_pts) || size(all_pts,1)<3
        return;
    end

    if strcmpi(roi_area_method,'convex')
        try
            k = convhull(all_pts(:,1), all_pts(:,2));
            A_px2 = polyarea(all_pts(k,1), all_pts(k,2));
            ROI_area_nm2 = A_px2 * (pixel_size_nm^2);
            ROI_area_um2 = ROI_area_nm2 / 1e6;
            ROI_regions = struct('X', all_pts(k,1), 'Y', all_pts(k,2));
            ROI_sampler.type  = 'convex';
            ROI_sampler.bbox  = [min(all_pts(:,1)) max(all_pts(:,1)) min(all_pts(:,2)) max(all_pts(:,2))];
            ROI_sampler.scale = NaN;
            return;
        catch
        end
    end

    % Rasterización multi-región
    xmin = floor(min(all_pts(:,1)) - 2);
    xmax = ceil(max(all_pts(:,1)) + 2);
    ymin = floor(min(all_pts(:,2)) - 2);
    ymax = ceil(max(all_pts(:,2)) + 2);

    s = max(1, round(mask_scale));
    W = max(1, round((xmax - xmin)*s) + 1);
    H = max(1, round((ymax - ymin)*s) + 1);

    jj = 1 + round((all_pts(:,1) - xmin) * s);
    ii = 1 + round((all_pts(:,2) - ymin) * s);
    jj = max(1, min(W, jj));  ii = max(1, min(H, ii));

    mask = false(H, W);
    mask(sub2ind([H W], ii, jj)) = true;

    rDil = max(0, round(mask_dilate_px * s));
    rCls = max(0, round(mask_close_px  * s));
    if rDil > 0, mask = imdilate(mask, strel('disk', rDil, 0)); end
    if rCls > 0, mask = imclose(mask, strel('disk', rCls, 0));  end

    pix_nm2 = (pixel_size_nm^2) / (s^2);
    min_pix = ceil( (min_region_area_um2 * 1e6) / pix_nm2 );
    if min_pix > 0, mask = bwareaopen(mask, min_pix, 8); end

    CC = bwconncomp(mask, 8);
    if CC.NumObjects > 0
        stats = regionprops(CC, 'Area');
        ROI_area_nm2 = sum([stats.Area]) * pix_nm2;
        ROI_area_um2 = ROI_area_nm2 / 1e6;

        B = bwboundaries(mask, 'noholes');
        for b = 1:numel(B)
            bc = B{b};
            if size(bc,1) < 3, continue; end
            ii_b = bc(:,1); jj_b = bc(:,2);
            xv = (jj_b - 1)/s + xmin;
            yv = (ii_b - 1)/s + ymin;
            ROI_regions(end+1).X = xv; %#ok<AGROW>
            ROI_regions(end).Y    = yv;
        end
        ROI_sampler.type  = 'mask-multi';
        ROI_sampler.bbox  = [xmin xmax ymin ymax];
        ROI_sampler.scale = s;
    end
end

function rawXY = read_float_pairs_xy(fpath)
% Interpreta el archivo como pares float32 [x y] consecutivos (fallback ROI).
    fid = fopen(fpath,'rb'); assert(fid>0, 'No se puede abrir el archivo.');
    data = fread(fid, inf, 'float32=>single'); fclose(fid);
    if mod(numel(data),2) ~= 0
        error('La longitud no es par; no se pueden formar pares [x y].');
    end
    rawXY = reshape(data, 2, []).';
end

function [mode_val, mode_count] = modal_bin_center(x, binw)
    x = x(isfinite(x));
    if isempty(x) || binw<=0
        mode_val = NaN; mode_count = 0; return;
    end
    xmin = min(x); xmax = max(x);
    if xmin==xmax
        mode_val = xmin; mode_count = numel(x); return;
    end
    start_edge = floor(xmin/binw)*binw;
    end_edge   = ceil(xmax/binw)*binw;
    edges = start_edge:binw:end_edge;
    [N,~] = histcounts(x, edges);
    if isempty(N) || all(N==0)
        mode_val = NaN; mode_count = 0; return;
    end
    [mode_count, idx] = max(N);
    mode_val = edges(idx) + binw/2;
end
