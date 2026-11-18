% MATLAB ROCK PHYSICS - INTERACTIVE COMPARISON WORKBENCH (test3.m)
% Run this standalone script in MATLAB with CoolProp configured.
% Original Fig1(a-g) and Fig2(a-d) are retained. Fig1(e): Sg=100%.
% Fig3/4: A, moving centers at fixed nuisance widths (phi/Sg and P/T).
% Fig5/6: B, changing nuisance widths at fixed centers (phi/Sg and P/T).
% Each NEW figure compares: 95% ellipse area, LDA Gaussian, 2D Gaussian,
% and LDA KDE overlap. No observation/model noise is added.
% All windows remain open. NO output directories or figure files are made.
% Workspace: cfg (settings), A/M (old data), maps/checks, figure_handles.
% For mathematics-only checks without fonts/CoolProp/figures:
%   setenv('RPT_TEST3_MODE','selftest'); test3;
%   setenv('RPT_TEST3_MODE','');       % restore ordinary interactive mode
clear; clc;

%% ---------- User-adjustable settings ------------------------------------
cfg = struct();
cfg.languages = {'CN', 'EN'};    % use {'CN'} or {'EN'} for fewer windows
cfg.compute_new_maps = true;
cfg.run_self_checks = true;
cfg.check_resolution = true;    % sample/method/grid probes; console + maps
cfg.check_direct_lut = true;    % compare representative LUT vs CoolProp
cfg.seed = 2026;
cfg.n_map = 1024;               % per gas/per pixel; same parameter samples
cfg.n_reference = 4096;         % Fig2 and the reference markers
cfg.batch_pixels = 16;          % vectorized forward-model block size
cfg.figure_position = [40 40 1280 850];

% Physical parameter support and reference settings.
cfg.phi_bounds = [0.05 0.30];
cfg.sg_bounds = [0.05 1.00];
cfg.T_bounds = [40 120];
cfg.P_bounds = [1 35]*1e6;
cfg.T_ref = 60; cfg.P_ref = 35e6;
cfg.phi_ref = 0.15; cfg.sg_ref = 0.60;
cfg.phi_fixed = 0.10; cfg.sg_fixed = 1.00; % Fig2/new maps; Fig1(e) stays Sg=100%
cfg.T_center = 80; cfg.P_center = 18e6;
cfg.sigma_phi = 0.025; cfg.sigma_sg = 0.120;
cfg.delta_T = 5; cfg.delta_P = 2e6; % A: adjustable UNIFORM half widths
% 'clipped' preserves the original min(max(normal,...),...) distribution.
% 'truncated' is an explicit alternative, applied to Fig2 AND every new map.
% Both use PRE-clipping/truncation Gaussian sigma, NOT the realized sigma.
cfg.normal_sampling = 'clipped';

% Actual calculation grids; fractions .00125/.005 = .125/.5 percentage points.
cfg.A_phi = 0.10:0.00125:0.25;
cfg.A_sg = 0.29:0.005:0.76;
cfg.A_T = 45:0.5:115;
cfg.A_P = (3:0.25:33)*1e6;
% Extra near-zero WIDTH coordinates resolve thin-cloud transition bands.
% These are real model evaluations, not visual interpolation. Rounding
% prevents effectively duplicate linear/logarithmic endpoints.
cfg.B_sigma_phi = unique(round(1e12*[0:0.0005:0.05, ...
    logspace(-6,log10(.005),32)])/1e12);
cfg.B_sigma_sg = unique(round(1e12*[0:0.0025:0.25, ...
    logspace(-5,log10(.025),32)])/1e12);
cfg.B_delta_T = 0:0.5:40;
cfg.B_delta_P = (0:0.25:17)*1e6;
cfg.contrast_phi = 0.05:0.00125:0.30;
cfg.contrast_sg = 0:0.005:1;
cfg.contrast_T = 40:0.5:120;
cfg.contrast_P = (1:0.25:35)*1e6;
cfg.lut_T = 35:1:125;
cfg.lut_P = (0.75:0.25:46)*1e6;

cfg.metric = struct('ellipse_points', 128, 'quadrature_order', 24, ...
    'kde_points', 256, 'compute_kde', true, 'rank_tol', 128*eps);
cfg.convergence_tolerance = 0.01; % absolute score: one percentage point
if cfg.run_self_checks || strcmp(getenv('RPT_TEST3_MODE'), 'selftest')
    run_test3_self_checks();
end
if strcmp(getenv('RPT_TEST3_MODE'), 'selftest'), return; end
validate_test3_settings(cfg);
check_figure_fonts();
rng(cfg.seed, 'twister');

%% ---------- 1. 矿物参数（石英与黏土，Mavko et al.）-----------------------
K_qtz    = 36.6e9;   % 石英体积模量 [Pa]
G_qtz    = 45.0e9;   % 石英剪切模量 [Pa]
rho_qtz  = 2650;     % 石英密度 [kg/m^3]
K_clay   = 21.0e9;   % 黏土体积模量 [Pa]
G_clay   = 7.0e9;    % 黏土剪切模量 [Pa]
rho_clay = 2600;     % 黏土密度 [kg/m^3]

%% ---------- 2. 固定岩石参数与计算范围设定 --------------------------------
alpha_sand = 0.12;   % 砂孔隙纵横比
alpha_clay = 0.035;  % 黏土孔隙纵横比
Vsh        = 0.30;   % 泥质体积含量
phi_fixed  = cfg.phi_fixed;   % 固定孔隙度 10%
salinity   = 0.035;  % 地层盐水矿化度 (35,000 ppm)

% 颜色样式规范
col_h2  = [0.52, 0.27, 0.75];  % 氢气 H2：紫色
col_ch4 = [0.82, 0.49, 0.08];  % 甲烷 CH4：金黄色
col_sw  = [0.11, 0.42, 0.75];  % 盐水 Brine：蓝色

%% ---------- 3. 固体基质模量 (VRH 平均) -----------------------------------
[K_solid, G_solid, rho_solid] = vrh([K_qtz, K_clay], [G_qtz, G_clay], ...
                                    [rho_qtz, rho_clay], [1 - Vsh, Vsh]);

%% ---------- 4. 干岩石骨架模量预积分 (DEM, Keys & Xu 2002) ----------------
N_dem = 1500;
phi_dem_axis = linspace(0, 0.35, N_dem);
dphi = phi_dem_axis(2) - phi_dem_axis(1);
K_dem = zeros(1, N_dem); G_dem = zeros(1, N_dem);
K_dem(1) = K_solid; G_dem(1) = G_solid;
K_curr = K_solid; G_curr = G_solid;
for k = 1:(N_dem - 1)
    p = phi_dem_axis(k);
    [Ps, Qs] = berryman_pq(K_curr, G_curr, 0, 0, alpha_sand);
    [Pc, Qc] = berryman_pq(K_curr, G_curr, 0, 0, alpha_clay);
    Peff = (1 - Vsh) * Ps + Vsh * Pc;
    Qeff = (1 - Vsh) * Qs + Vsh * Qc;
    
    K_curr = K_curr - K_curr * Peff * dphi / (1 - p);
    G_curr = G_curr - G_curr * Qeff * dphi / (1 - p);
    K_dem(k + 1) = K_curr;
    G_dem(k + 1) = G_curr;
end
Kd_10 = interp1(phi_dem_axis, K_dem, phi_fixed, 'pchip');
Gd_10 = interp1(phi_dem_axis, G_dem, phi_fixed, 'pchip');

%% ========================================================================
%  画布 1：理想岩石物理模型下的 H2-CH4 弹性响应差异
%  不引入观测误差或误差协方差；仅比较正演属性的相对变化。
%  对两种量纲不同的属性分别计算百分比差异，再做等权 RMS 组合：
%    C_Ip = 100*abs(Ip_H2-Ip_CH4) / (0.5*(Ip_H2+Ip_CH4))
%    C_R  = 100*abs(R_H2-R_CH4)   / (0.5*(R_H2+R_CH4)), R=Vp/Vs
%    C_joint = sqrt((C_Ip^2+C_R^2)/2)         [单位：%]
%  C_joint 是无观测误差的物理响应差异指标，不是 Mahalanobis 距离，
%  不能直接换算为检测概率、误判率或统计显著性。
%  ========================================================================
T_rpt = cfg.T_ref;         % 温度基准 [degC]
P_rpt = cfg.P_ref;       % 孔隙压力基准 [Pa]

phi_vec = cfg.contrast_phi; Sg_vec = cfg.contrast_sg;
N_phi = numel(phi_vec); N_sg = numel(Sg_vec);
[PHI_grid, SG_grid] = meshgrid(phi_vec, Sg_vec);
Kd_grid = interp1(phi_dem_axis, K_dem, PHI_grid, 'pchip');
Gd_grid = interp1(phi_dem_axis, G_dem, PHI_grid, 'pchip');

[Ip_c_grid, VpVs_c_grid] = calc_rpt_phi_sg(Kd_grid, Gd_grid, PHI_grid, SG_grid, ...
    T_rpt, P_rpt, 'Methane', salinity, K_solid, rho_solid);
[Ip_h_grid, VpVs_h_grid] = calc_rpt_phi_sg(Kd_grid, Gd_grid, PHI_grid, SG_grid, ...
    T_rpt, P_rpt, 'Hydrogen', salinity, K_solid, rho_solid);
[C_joint_phi, C_ip_phi, C_ratio_phi] = calc_ideal_contrast( ...
    Ip_h_grid, VpVs_h_grid, Ip_c_grid, VpVs_c_grid);

% CoolProp 查找表保留，用于大规模 P-T 组合与画布 2 的随机采样。
% LUT 的边界覆盖当前 Monte Carlo 压力、温度设置下的大部分样本；
% 超出插值范围的样本会显式报错，不作静默外推。
fprintf('Building refined CoolProp LUT (%d x %d states per gas)...\n', ...
    numel(cfg.lut_T), numel(cfg.lut_P));
gas_lut = build_coolprop_lut(cfg.lut_T, cfg.lut_P);

% P-T signal map: refined physical grid, fixed Sg=100%.
P_pt = cfg.contrast_P; T_pt = cfg.contrast_T;
[PP_pt, TT_pt] = meshgrid(P_pt, T_pt);
phi_pt_grid = phi_fixed * ones(size(PP_pt));
sg_pt_grid  = ones(size(PP_pt));
Kd_pt_grid  = Kd_10 * ones(size(PP_pt));
Gd_pt_grid  = Gd_10 * ones(size(PP_pt));
[Ip_c_pt_map, VpVs_c_pt_map] = calc_rpt_lut(Kd_pt_grid, Gd_pt_grid, ...
    phi_pt_grid, sg_pt_grid, TT_pt, PP_pt, 'Methane', salinity, K_solid, rho_solid, gas_lut);
[Ip_h_pt_map, VpVs_h_pt_map] = calc_rpt_lut(Kd_pt_grid, Gd_pt_grid, ...
    phi_pt_grid, sg_pt_grid, TT_pt, PP_pt, 'Hydrogen', salinity, K_solid, rho_solid, gas_lut);
C_pt = calc_ideal_contrast(Ip_h_pt_map, VpVs_h_pt_map, ...
    Ip_c_pt_map, VpVs_c_pt_map);

%% ========================================================================
%  Compute once: reuse identical numerical results in Chinese and English.
%  Reuse one physical computation for both languages; no export is performed.
%  Figure 1 panels: a=Ip, b=Vp/Vs, c=RPT, d=joint contrast,
%  e=P-T map, f=pressure sweep, g=temperature sweep.
%  Panels a,b,d,e,f,g: 2.90 x 3.60 in; panel c: 6.40 x 3.60 in.
%  ========================================================================
A = struct();
A.phi_vec = phi_vec; A.Sg_vec = Sg_vec;
A.PHI_grid = PHI_grid; A.SG_grid = SG_grid;
A.C_joint_phi = C_joint_phi; A.C_ip_phi = C_ip_phi;
A.C_ratio_phi = C_ratio_phi; A.C_pt = C_pt;
A.P_pt = P_pt; A.T_pt = T_pt;
A.PP_pt = PP_pt; A.TT_pt = TT_pt;
A.phi_dem_axis = phi_dem_axis; A.K_dem = K_dem; A.G_dem = G_dem;
A.Kd_10 = Kd_10; A.Gd_10 = Gd_10;
A.T_rpt = T_rpt; A.P_rpt = P_rpt;
A.salinity = salinity; A.K_solid = K_solid;
A.rho_solid = rho_solid; A.phi_fixed = phi_fixed;
A.col_h2 = col_h2; A.col_ch4 = col_ch4; A.col_sw = col_sw;
A.cfg = cfg;
% Evaluate RPT mesh and sensitivity curves only once for both languages.
A.V = precompute_visual_curves(A);

%% ========================================================================
%  Figure 2: paired parameter ensembles, not measurement noise.
%  This random sampling is executed ONCE so the two language versions match.
%  ========================================================================
N_mc = cfg.n_reference;
U_reference = shared_sample_design(N_mc, cfg.seed);
phi_dense_rpt = cfg.contrast_phi;
[Ip_sw_line, VpVs_sw_line] = calc_rpt_lut( ...
    interp1(phi_dem_axis, K_dem, phi_dense_rpt, 'pchip'), ...
    interp1(phi_dem_axis, G_dem, phi_dense_rpt, 'pchip'), ...
    phi_dense_rpt, zeros(size(phi_dense_rpt)), T_rpt, P_rpt, ...
    'Methane', salinity, K_solid, rho_solid, gas_lut);

% Row 1: vary ONLY porosity and gas saturation. T and P are STRICTLY fixed.
[phi_mc1, Sg_mc1] = nuisance_samples('B_phi', ...
    [cfg.sigma_phi cfg.sigma_sg], cfg, U_reference);
T_mc1   = T_rpt * ones(N_mc, 1);   % fixed 60 degC (NO stochastic perturbation)
P_mc1   = P_rpt * ones(N_mc, 1);   % fixed 35 MPa (NO stochastic perturbation)
assert(all(T_mc1 == T_rpt) && all(P_mc1 == P_rpt), ...
    'RPT:Row1TPNotFixed', 'Figure 2 top row must be exactly constant T and P.');
Kd_mc1  = interp1(phi_dem_axis, K_dem, phi_mc1, 'pchip');
Gd_mc1  = interp1(phi_dem_axis, G_dem, phi_mc1, 'pchip');
[Ip_c_r1, VpVs_c_r1] = calc_rpt_lut(Kd_mc1, Gd_mc1, phi_mc1, Sg_mc1, ...
    T_mc1, P_mc1, 'Methane', salinity, K_solid, rho_solid, gas_lut);
[Ip_h_r1, VpVs_h_r1] = calc_rpt_lut(Kd_mc1, Gd_mc1, phi_mc1, Sg_mc1, ...
    T_mc1, P_mc1, 'Hydrogen', salinity, K_solid, rho_solid, gas_lut);

% Row 2: fixed phi=10% and Sg=100%; ONLY T and P vary.
phi_mc2 = phi_fixed*ones(N_mc, 1);
cfg_reference_tp = reference_tp_settings(cfg);
[phi_mc2, Sg_mc2, T_mc2, P_mc2] = nuisance_samples('B_tp', ...
    [diff(cfg.P_bounds)/2 diff(cfg.T_bounds)/2], cfg_reference_tp, U_reference);
Kd_mc2  = Kd_10*ones(N_mc, 1);
Gd_mc2  = Gd_10*ones(N_mc, 1);
[Ip_c_r2, VpVs_c_r2] = calc_rpt_lut(Kd_mc2, Gd_mc2, phi_mc2, Sg_mc2, ...
    T_mc2, P_mc2, 'Methane', salinity, K_solid, rho_solid, gas_lut);
[Ip_h_r2, VpVs_h_r2] = calc_rpt_lut(Kd_mc2, Gd_mc2, phi_mc2, Sg_mc2, ...
    T_mc2, P_mc2, 'Hydrogen', salinity, K_solid, rho_solid, gas_lut);

M = struct();
M.Ip_c_r1 = Ip_c_r1; M.VpVs_c_r1 = VpVs_c_r1;
M.Ip_h_r1 = Ip_h_r1; M.VpVs_h_r1 = VpVs_h_r1;
M.Ip_c_r2 = Ip_c_r2; M.VpVs_c_r2 = VpVs_c_r2;
M.Ip_h_r2 = Ip_h_r2; M.VpVs_h_r2 = VpVs_h_r2;
M.Ip_sw_line = Ip_sw_line; M.VpVs_sw_line = VpVs_sw_line;
M.col_h2 = col_h2; M.col_ch4 = col_ch4; M.col_sw = col_sw;
M.cfg = cfg;
reference_metric = cfg.metric; reference_metric.compute_kde = true;
M.stats_r1 = ensemble_separability([Ip_c_r1 VpVs_c_r1], ...
    [Ip_h_r1 VpVs_h_r1], reference_metric);
M.stats_r2 = ensemble_separability([Ip_c_r2 VpVs_c_r2], ...
    [Ip_h_r2 VpVs_h_r2], reference_metric);
print_reference_scores('Fig2 phi/Sg', M.stats_r1);
print_reference_scores('Fig2 T/P', M.stats_r2);

%% ---------- Interactive windows: original panels remain available --------
% This script intentionally contains NO print/export/save/mkdir operations.
% Handles and numeric data stay in the base workspace for manual adjustments.
slide_height_in = 9;
panel_height_in = (2/5)*slide_height_in;
fig1_portrait_size_in = [2.90, panel_height_in];
landscape_size_in = [(16/9)*panel_height_in, panel_height_in];
lang_list = cfg.languages;
fig1_panels = 'abcdefg';
fig2_panels = 'abcd';
figure_handles = struct('Fig1', [], 'Fig2', [], 'comparisons', []);
figure_handles.Fig1 = gobjects(numel(lang_list), numel(fig1_panels));
figure_handles.Fig2 = gobjects(numel(lang_list), numel(fig2_panels));
figure_handles.comparisons = gobjects(numel(lang_list), 4);
for il = 1:numel(lang_list)
    lang = lang_list{il};
    for k = 1:numel(fig1_panels)
        fig = render_ideal_panel(A, lang, fig1_panels(k), ...
            fig1_portrait_size_in, landscape_size_in);
        lock_final_panel_typography(fig);
        figure_handles.Fig1(il,k) = fig;
        drawnow;
    end
    for k = 1:numel(fig2_panels)
        fig = render_mc_panel(M, lang, k, fig2_panels(k), landscape_size_in);
        lock_final_panel_typography(fig);
        figure_handles.Fig2(il,k) = fig;
        drawnow;
    end
end

%% ---------- Supplementary maps: A/B, all four overlap measures ------------
% S_ellipse is an AREA HEURISTIC, not a classification probability.
% Gaussian/KDE S=1-OVL=1-2Pe assumes equal class priors, for that fitted model.
% S>=.99 means near-complete separation; never round .99 up to exact1.
% Degenerate zero-area ellipses/undefined LDA directions are shown in gray.
maps = cell(1,4); checks = struct();
if cfg.compute_new_maps
    U_map = shared_sample_design(cfg.n_map, cfg.seed);
    if cfg.check_direct_lut
        checks.lut = check_lut_references(A, gas_lut, cfg, U_reference);
    end
    map_kinds = {'A_phi','A_tp','B_phi','B_tp'};
    for im = 1:numel(map_kinds)
        maps{im} = build_separability_map(map_kinds{im}, A, gas_lut, cfg, U_map);
        if cfg.check_resolution
            maps{im}.checks = check_separability_resolution(maps{im}, A, gas_lut, cfg);
        end
        for il = 1:numel(lang_list)
            figure_handles.comparisons(il,im) = ...
                render_separability_comparison(maps{im}, lang_list{il}, cfg);
            drawnow;
        end
    end
end
fprintf(['Interactive figures are ready; no figure files were written.\n' ...
    'Edit cfg at the top of test3.m; inspect maps, checks and figure_handles.\n']);

%% ========================================================================
%  STANDALONE FIGURE 1, one panel per figure. No tiledlayout or title().
%  Letter IDs are figure-space annotations; data axes stay unobstructed.
%  Scientific results A / A.V are calculated once for CN and EN.
%  ========================================================================
function fig = create_panel_figure(fig_name, page_size_in)
% True physical page dimensions; PDF and SVG inherit these inches.
% On-screen pixels are just a drawing resolution, NOT the final output size.
if nargin < 2 || isempty(page_size_in)
    error('RPT:MissingPanelSize', ...
        'Every standalone panel must specify its physical page size.');
end
assert(isnumeric(page_size_in) && numel(page_size_in)==2 && ...
    all(isfinite(page_size_in)) && all(page_size_in > 0), ...
    'RPT:InvalidPaperSize', 'Paper size must be [width height] in inches.');
if page_size_in(1) < 4
    screen_ppi = 230;  % comfortably readable window for the small panels
else
    screen_ppi = 120;  % original 1920 x 1080 screen window
end
size_px = round(screen_ppi * page_size_in);
fig = figure('Name', fig_name, 'NumberTitle', 'off', 'Visible', 'on', ...
    'Color', 'w', 'Units', 'pixels', ...
    'Position', [40 40 size_px(1) size_px(2)], 'Renderer', 'painters');
setappdata(fig, 'PhysicalExportSizeInches', page_size_in);
end

function fig = render_ideal_panel(A, lang, id, portrait_size_in, landscape_size_in)
compact = any(id == 'abdefg');
if compact
    fig = create_panel_figure(sprintf('Figure 1(%s) - %s', id, lang), ...
        portrait_size_in);
else
    fig = create_panel_figure(sprintf('Figure 1(%s) - %s', id, lang), ...
        landscape_size_in);
end
ax = axes('Parent', fig, 'Units', 'normalized', ...
    'Position', [.135 .17 .715 .675], 'Color', 'white');
V = A.V;
cmap = turbo(256);
common_max = max([A.C_joint_phi(:); A.C_ip_phi(:); A.C_ratio_phi(:)]);
common_max = max(common_max, 1e-9);

switch id
    case 'a'
        plot_contrast_map(ax, 100*A.phi_vec, 100*A.Sg_vec, ...
            100*A.PHI_grid, 100*A.SG_grid, A.C_ip_phi, ...
            [0 common_max], cmap, lang, ...
            'Iₚ 相对差异 (%)', 'Iₚ relative contrast (%)', true);
    case 'b'
        plot_contrast_map(ax, 100*A.phi_vec, 100*A.Sg_vec, ...
            100*A.PHI_grid, 100*A.SG_grid, A.C_ratio_phi, ...
            [0 common_max], cmap, lang, ...
            'Vₚ/Vₛ 相对差异 (%)', 'Vₚ/Vₛ relative contrast (%)', true);
    case 'c'
        % Single 6.40 x 3.60 in RPT panel, without a zoom inset.
        % Use the full central plotting region; reserve a bottom band for the key.
        set(ax, 'Position', [.145 .170 .785 .685]);
        hold(ax, 'on');
        draw_rpt_mesh(ax, V, A.col_ch4, A.col_h2, A.col_sw);
        % Reference-point markers and connecting line removed.
        style_axes(ax, 18);
        apply_axis_labels(ax, lang, ...
            '纵波阻抗 Iₚ (MRayl)', 'P-wave impedance Iₚ (MRayl)', ...
            '纵横波速比 Vₚ/Vₛ', 'Vₚ/Vₛ');
        grid(ax, 'on');
        annotate_rpt_grid(ax, V, A.col_sw, A.col_ch4);
        draw_rpt_key_inside(ax, lang, V, A.col_sw, A.col_ch4, A.col_h2);
        format_landscape_rpt_axes(fig, ax);
        % No secondary axes / local-detail inset in Figure 1(c).
    case 'd'
        plot_contrast_map(ax, 100*A.phi_vec, 100*A.Sg_vec, ...
            100*A.PHI_grid, 100*A.SG_grid, A.C_joint_phi, ...
            [0 common_max], cmap, lang, ...
            'RMS 相对差异 (%)', 'RMS relative contrast (%)', true);
        add_portrait_condition_header(fig, sprintf('T = %g °C, P = %g MPa', A.T_rpt, A.P_rpt/1e6));
    case 'e'
        hold(ax, 'on');
        draw_vector_heatmap(ax, A.PP_pt/1e6, A.TT_pt, A.C_pt);
        set(ax, 'YDir', 'normal');
        colormap(ax, cmap);
        pt_min = min(A.C_pt(:)); pt_max = max(A.C_pt(:));
        clim(ax, [0 max(1e-9, pt_max*1.035)]);
        if pt_max > pt_min
            levels = linspace(pt_min, pt_max, 5);
            overlay_value_contours(ax, A.PP_pt/1e6, A.TT_pt, ...
                A.C_pt, levels(2:4), true);
        end
        style_axes(ax, 17.3);
        apply_axis_labels(ax, lang, ...
            '孔隙压力 P (MPa)', 'Pore pressure P (MPa)', ...
            '温度 T (°C)', 'Temperature T (°C)');
        xlim(ax, [min(A.P_pt) max(A.P_pt)]/1e6); ylim(ax, [min(A.T_pt) max(A.T_pt)]);
        style_colorbar(ax, lang, 'RMS 相对差异 (%)', ...
            'RMS relative contrast (%)');
        add_portrait_condition_header(fig, sprintf('φ = %g%%, S_g = 100%%', 100*A.phi_fixed));
    case {'f', 'g'}
        plot_ideal_sweep(ax, fig, A, lang, id);
    otherwise
        error('RPT:InvalidFigure1Panel', 'Unknown Figure 1 panel: %s.', id);
end
if compact
    % Make the ACTUAL data axes vertical, not just the white page.
    % Reposition the colorbar and resize fonts for final 3.6-inch PPT use.
    format_portrait_ideal_panel(fig, ax, id);
end
add_panel_id(fig, ['(' id ')']);
if id == 'c'
    assert(numel(findall(fig, 'Type', 'Axes')) == 1, ...
        'RPT:UnexpectedInset', 'Figure 1(c) must contain only one main axes.');
end
assert_no_panel_titles(fig);
end

function plot_ideal_sweep(ax, fig, A, lang, id)
V = A.V;
colors = [.06 .53 .44; .15 .41 .77; .83 .27 .30];
styles = {'--', '-.', '-'};
peak = 0;
hold(ax, 'on');
for i = 1:numel(V.Sg_list)
    if id == 'f'
        xx = V.P_sweep/1e6;
        yy = V.C_p{i};
    else
        xx = V.T_sweep;
        yy = V.C_t{i};
    end
    plot(ax, xx, yy, styles{i}, ...
        'Color', colors(i,:), 'LineWidth', 2.6, ...
        'HandleVisibility', 'off');
    peak = max(peak, max(yy(:)));
end
style_axes(ax, 17.5);
if id == 'f'
    apply_axis_labels(ax, lang, ...
        '孔隙压力 P (MPa)', 'Pore pressure P (MPa)', ...
        'RMS 相对差异 (%)', 'RMS relative contrast (%)');
    xlim(ax, [min(V.P_sweep) max(V.P_sweep)]/1e6);
    add_small_condition_note(ax, {sprintf('φ = %g%%',100*A.phi_fixed), ...
        sprintf('T = %g °C',A.T_rpt)}, 'southeast');
else
    apply_axis_labels(ax, lang, ...
        '温度 T (°C)', 'Temperature T (°C)', ...
        'RMS 相对差异 (%)', 'RMS relative contrast (%)');
    xlim(ax, [min(V.T_sweep) max(V.T_sweep)]);
    add_small_condition_note(ax, {sprintf('φ = %g%%',100*A.phi_fixed), ...
        sprintf('P = %g MPa',A.P_rpt/1e6)}, 'southeast');
end
if peak > 0, ylim(ax, [0 1.12*peak]); end
grid(ax, 'on');
% Move the S_g key inside the axes so the narrow portrait-page header stays
% clean and the legend cannot collide with the panel identifier.
draw_sweep_key_inside(ax, colors, styles);
end

function draw_sweep_key_inside(ax, colors, styles)
% Figure 1(f,g): compact in-axes legend in the upper-left blank zone.
levels = [20 50 100];
xl = xlim(ax); yl = ylim(ax);
dx = diff(xl); dy = diff(yl);
x1 = xl(1) + 0.08*dx;
x2 = xl(1) + 0.28*dx;
text_x = xl(1) + 0.31*dx;
yf = [0.92 0.84 0.76];
for k = 1:numel(levels)
    yy = yl(1) + yf(k)*dy;
    plot(ax, [x1 x2], [yy yy], styles{k}, 'Color', colors(k,:), ...
        'LineWidth', 1.85, 'Clipping', 'on', 'HandleVisibility', 'off');
    [s, intrp] = scientific_label(sprintf('S_g = %d%%', levels(k)));
    text(ax, text_x, yy, s, ...
        'Interpreter', intrp, 'FontName', 'Times New Roman', ...
        'FontSize', 10.4, 'FontWeight', 'bold', ...
        'Color', [.17 .20 .24], 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle', 'BackgroundColor', 'w', ...
        'Margin', 0.25, 'Clipping', 'on');
end
end

function add_portrait_condition_header(fig, condition)
% State descriptions are not subplot titles: keep them in the slim header,
% outside the heatmap axes, so no numerical contour is covered.
[s, intrp] = scientific_label(condition);
annotation(fig, 'textbox', [.235 .873 .710 .047], ...
    'String', s, 'Interpreter', intrp, ...
    'FontName', 'Times New Roman', ...
    'FontSize', portrait_header_font_size(intrp), ...
    'FontWeight', 'bold', 'Color', [.17 .21 .25], ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'middle', 'FitBoxToText', 'off', 'Margin', 0);
end

function fs = portrait_header_font_size(interpreter)
% The parameter header in (d) is plain text, whereas (e) contains S_{g}
% and requires MATLAB TeX. Compensate for the TeX/PDF glyph-size reduction
% so both headers appear at the same PHYSICAL size after PDF/SVG export.
fs = 10.5;
if strcmpi(char(interpreter), 'tex')
    fs = fs * (11.4 / 9.0);
end
end

function format_portrait_ideal_panel(fig, ax, id)
% Paper size: 2.90 x 3.60 in.  Every data axes must also be portrait.
% Layout is in normalized page units; positions are measured after drawing
% to avoid colorbar automatic layout taking away space unexpectedly.
page_in = getappdata(fig, 'PhysicalExportSizeInches');
assert(abs(page_in(2)-3.6) < 1e-8 && page_in(1)<page_in(2), ...
    'RPT:InvalidPortraitHeight', 'Expected portrait page height = 3.60 inches.');
set(ax, 'Units', 'normalized');
if any(id == 'abde')
    % Heatmap: reserve the right-hand column for a narrow colorbar + label.
    % Axis width = 0.480*2.90 = 1.392 in;
    % axis height = 0.700*3.60 = 2.520 in.
    ax_position = [.235 .150 .480 .700];
    cb_position = [.766 .150 .038 .700];
else
    % Sensitivity lines: keep exactly the SAME main-axes rectangle as
    % (a,b,d,e), even though these panels do not require a colorbar.
    % This ensures all six portrait data plots align when placed in PPT.
    % Axis width = 0.480*2.90 = 1.392 in;
    % axis height = 0.700*3.60 = 2.520 in.
    ax_position = [.235 .150 .480 .700];
    cb_position = [];
end
set(ax, 'Position', ax_position);
actual_w_in = ax_position(3)*page_in(1);
actual_h_in = ax_position(4)*page_in(2);
assert(actual_h_in > actual_w_in, 'RPT:AxisNotPortrait', ...
    'Figure 1(%s) data axes must be taller than wide.', id);

% Native-size typography: tick labels 10.3 pt, X/Y labels visually
% 11.5 pt in PDF/SVG even for TeX-rendered scientific subscripts.
set(ax, 'FontSize', 10.3, 'LineWidth', 1.15, ...
    'TickLength', [.015 .018]);
texts = findall(ax, 'Type', 'Text');
for k = 1:numel(texts)
    set(texts(k), 'FontSize', 10.4);
end
set(get(ax, 'XLabel'), 'FontUnits', 'points', ...
        'FontSize', axis_label_font_size(get(get(ax, 'XLabel'), 'Interpreter')), ...
        'FontWeight', 'bold');
set(get(ax, 'YLabel'), 'FontUnits', 'points', ...
        'FontSize', axis_label_font_size(get(get(ax, 'YLabel'), 'Interpreter')), ...
        'FontWeight', 'bold');

if any(id == 'abde')
    cb = findall(fig, 'Type', 'ColorBar');
    assert(numel(cb) == 1, 'RPT:MissingColorbar', ...
        'Figure 1(%s) requires exactly one colorbar.', id);
    set(cb, 'Units', 'normalized', 'Position', cb_position, ...
        'FontSize', 10.3, 'LineWidth', 0.95);
    set(cb.Label, 'FontUnits', 'points', ...
        'FontSize', colorbar_label_font_size(get(cb.Label, 'Interpreter')), ...
        'FontWeight', 'bold');
    if id == 'e'
        set(ax,'XTickMode','auto','YTickMode','auto');
    else
        set(ax,'XTickMode','auto','YTickMode','auto');
    end
elseif id == 'f'
    set(ax,'XTickMode','auto');
elseif id == 'g'
    set(ax,'XTickMode','auto');
end
% Ensure neither automatic axes nor colorbar layout changes the intended
% data-axes geometry at the end of styling.
set(ax, 'Position', ax_position);
end

function format_landscape_rpt_axes(fig, ax)
% Figure 1(c): one RPT axes; use the EXACT same font hierarchy as
% Figure 1(a,b,d,e,f,g): ticks 10.3, labels 11.5, annotations 10.4 pt.
page = getappdata(fig, 'PhysicalExportSizeInches');
assert(is_landscape_panel_size(page), 'RPT:InvalidRPTLandscape', ...
    'Figure 1(c) must be 6.40 x 3.60 inches.');
set(ax, 'FontSize', 10.3, 'LineWidth', 1.15, ...
    'TickLength', [.015 .018]);
texts = findall(ax, 'Type', 'Text');
for k = 1:numel(texts)
    set(texts(k), 'FontSize', 10.4);
end
set(get(ax, 'XLabel'), 'FontUnits', 'points', ...
        'FontSize', axis_label_font_size(get(get(ax, 'XLabel'), 'Interpreter')), ...
        'FontWeight', 'bold');
set(get(ax, 'YLabel'), 'FontUnits', 'points', ...
        'FontSize', axis_label_font_size(get(get(ax, 'YLabel'), 'Interpreter')), ...
        'FontWeight', 'bold');
end

function tf = is_landscape_panel_size(page)
tf = isnumeric(page) && numel(page)==2 && all(isfinite(page)) && ...
    abs(page(1)-6.4)<1e-8 && abs(page(2)-3.6)<1e-8 && ...
    page(1)>page(2);
end

function plot_contrast_map(ax, xv, yv, X, Y, C, clim_bounds, ...
        cmap, lang, cn_colorbar, en_colorbar, label_contours)
hold(ax, 'on');
draw_vector_heatmap(ax, X, Y, C);
set(ax, 'YDir', 'normal');
colormap(ax, cmap); clim(ax, clim_bounds);
level_candidates = [.25 .5 1 2 4 6 8 10 15 20 30 40];
levels = level_candidates(level_candidates > min(C(:)) & ...
                           level_candidates < max(C(:)));
if ~isempty(levels)
    if numel(levels) == 1, levels = [levels levels]; end
    [cc, hh] = contour(ax, X, Y, C, levels, '-', ...
        'Color', [.99 .99 .99], 'LineWidth', 1.45);
    if label_contours
        apply_contour_labels(cc, hh);
    end
end
style_axes(ax, 17.0);
apply_axis_labels(ax, lang, ...
    '孔隙度 φ (%)', 'Porosity φ (%)', ...
    '气体饱和度 S_g (%)', 'Gas saturation S_g (%)');
xlim(ax, [min(xv) max(xv)]); ylim(ax, [min(yv) max(yv)]);
style_colorbar(ax, lang, cn_colorbar, en_colorbar);
end

function overlay_value_contours(ax, X, Y, C, levels, label_contours)
% White contour lines with optional numeric labels, matching panel (d).
if nargin < 6, label_contours = true; end
levels = levels(:)';
levels = levels(isfinite(levels));
if isempty(levels), return; end
if numel(levels) == 1, levels = [levels levels]; end
[cc, hh] = contour(ax, X, Y, C, levels, '-', ...
    'Color', [.99 .99 .99], 'LineWidth', 1.35);
if label_contours
    apply_contour_labels(cc, hh, '%.2f');
end
end

function apply_contour_labels(cc, hh, number_format)
% Clearly formatted labels; e.g. show 0.10 rather than
% 0.098888 on the P/T diagram, without rounding the true contour curves.
if nargin < 3, number_format = ''; end
% Pass the FINAL font size to clabel so it reserves the correct gap
% within the white contour. Enlarging the text only afterwards would
% leave the gap too short and cause white lines to cross the digits.
% Tight portrait panels need physically shorter contour-label spacing;
% otherwise the original 485-pt interval can suppress all numeric labels.
page_in = getappdata(ancestor(hh, 'figure'), 'PhysicalExportSizeInches');
if ~isempty(page_in) && page_in(1) < page_in(2)
    spacing_pt = 135;
    label_font_pt = 10.4;
else
    spacing_pt = 485;
    label_font_pt = 15.2;
end
% Automatic labels managed by a Contour object do NOT support text-handle
% output from clabel(C,h,...). Requesting it creates a MATLAB warning and
% returns no usable handles. Keep the contour object to preserve rotated
% labels, original line gaps, and the intended label spacing, but do not
% request an output argument.
%
% MATLAB R2022b+ supports Contour.LabelFormat. Set it on the contour
% object instead of trying to modify inaccessible managed text handles.
% Older versions will still draw the labels with default numeric precision.
if ~isempty(number_format) && isprop(hh, 'LabelFormat')
    set(hh, 'LabelFormat', number_format);
end
clabel(cc, hh, 'LabelSpacing', spacing_pt, ...
    'FontName', 'Times New Roman', 'FontSize', label_font_pt, ...
    'FontWeight', 'bold', 'Color', [1 1 1], ...
    'Interpreter', 'none');
end

function annotate_rpt_grid(ax, V, col_brine, col_saturation)
% Numerical annotations reference calculated RPT nodes, never estimates.
% Porosity: position on brine-only reference curve (S_g = 0).
% Gas saturation: points at phi = 30% on CH4 constant-S_g curves.
% All numeric text is TNR, independent of the figure language.
old_x = xlim(ax); old_y = ylim(ax);
wx = diff(old_x); wy = diff(old_y);
xlim(ax, [old_x(1)-.22*wx, old_x(2)+.090*wx]);
ylim(ax, [old_y(1)-.11*wy, old_y(2)+.095*wy]);
full_x = xlim(ax); full_y = ylim(ax);
wx = diff(full_x); wy = diff(full_y);

% Selected porosity nodes: reference brine curve is common to H2/CH4.
B = V.Sat_C{1};
phi_to_label = V.phi_nodes([1 2 4 6]); % label 5, 10, 20, 30%; all grid curves retained
for k = 1:numel(phi_to_label)
    value = phi_to_label(k);
    px = interp1(V.phi_dense, B(:,1), value, 'linear');
    py = interp1(V.phi_dense, B(:,2), value, 'linear');
    plot(ax, px, py, 'o', 'MarkerSize', 4.1, ...
        'MarkerFaceColor', col_brine, 'MarkerEdgeColor', 'w', ...
        'LineWidth', .6, 'HandleVisibility', 'off');
    % Move labels into a clear zone above/right of the brine branch.
    % Each node has an independent, deliberately offset callout; this
    % prevents labels from sitting on the blue curve or crossing each other.
    % Labels are ordered phi = [5, 10, 20, 30] percent.
    label_dx = [.015, -.045, .055, .030] * wx;
    label_dy = [.095, .105, .115, .125] * wy;
    tx = px + label_dx(k);
    ty = py + label_dy(k);
    tx = min(max(tx, full_x(1)+.075*wx), full_x(2)-.095*wx);
    ty = min(ty, full_y(2)-.065*wy);
    % Thin blue leader lines anchor each shifted label to its true node.
    plot(ax, [px tx], [py ty-.018*wy], '-', ...
        'Color', .72*col_brine+.28*[1 1 1], 'LineWidth', 1.0, ...
        'HandleVisibility', 'off');
    text(ax, tx, ty, sprintf('φ=%d%%', round(100*value)), ...
        'Interpreter', 'none', 'FontName', 'Times New Roman', ...
        'FontSize', 10.4, 'FontWeight', 'bold', ...
        'Color', [.08 .32 .61], 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'BackgroundColor', 'w', 'Margin', 1.0, 'Clipping', 'on');
end

% Four gas-saturation annotations on CH4 contours, phi=30% endpoints.
sel = 2:numel(V.Sg_nodes); % 25%, 50%, 75%, 100%; S_g = 0 is brine curve
p = zeros(numel(sel),2);
for k = 1:numel(sel)
    C = V.Sat_C{sel(k)};
    p(k,:) = C(end,:);
end
label_y = p(:,2);
[~, order] = sort(label_y);
min_gap = .072*wy;
for k = 2:numel(order)
    curr = order(k); prev = order(k-1);
    label_y(curr) = max(label_y(curr), label_y(prev)+min_gap);
end
% Keep the last leader and all text inside plotting limits.
over = max(label_y)-(full_y(2)-.055*wy);
if over > 0, label_y = label_y-over; end
under = (full_y(1)+.075*wy)-min(label_y);
if under > 0, label_y = label_y+under; end
label_x = max(full_x(1)+.115*wx, min(p(:,1))-.12*wx);
for k = 1:numel(sel)
    plot(ax, p(k,1), p(k,2), 's', ...
        'MarkerSize', 4.2, 'MarkerFaceColor', col_saturation, ...
        'MarkerEdgeColor', 'w', 'LineWidth', .6, ...
        'HandleVisibility', 'off');
    plot(ax, [label_x+.017*wx p(k,1)], [label_y(k) p(k,2)], '-', ...
        'Color', .75*col_saturation+.25*[1 1 1], ...
        'LineWidth', 1.2, 'HandleVisibility', 'off');
    [sg_text, sg_intrp] = scientific_label( ...
        sprintf('S_g = %d%%', round(100*V.Sg_nodes(sel(k)))));
    text(ax, label_x, label_y(k), sg_text, ...
        'FontName', 'Times New Roman', 'Interpreter', sg_intrp, ...
        'FontSize', 16.2, 'FontWeight', 'bold', ...
        'Color', [.52 .31 .06], 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'middle', 'BackgroundColor', 'w', ...
        'Margin', .5, 'Clipping', 'on');
end
end

function fig = render_mc_panel(M, lang, panel_index, id, landscape_size_in)
% a/c: scatter + marginal KDE; b/d: LDA projection with KDE overlap.
% T/P and phi/Sg distributions are unchanged from the numerical model.
fig = create_panel_figure(sprintf('Figure 2(%s) - %s', id, lang), ...
    landscape_size_in);
if panel_index <= 2
    Ip_c = M.Ip_c_r1; R_c = M.VpVs_c_r1;
    Ip_h = M.Ip_h_r1; R_h = M.VpVs_h_r1;
    brine = true;
    stat = M.stats_r1;
    [cn_note,en_note] = fig2_condition_text(M.cfg,1);
else
    Ip_c = M.Ip_c_r2; R_c = M.VpVs_c_r2;
    Ip_h = M.Ip_h_r2; R_h = M.VpVs_h_r2;
    brine = false;
    stat = M.stats_r2;
    [cn_note,en_note] = fig2_condition_text(M.cfg,2);
end
Ip_c = Ip_c(:); R_c = R_c(:);
Ip_h = Ip_h(:); R_h = R_h(:);
if mod(panel_index, 2) == 1
    draw_mc_scatter_panel(fig, M, lang, ...
        Ip_c, R_c, Ip_h, R_h, brine);
    add_mc_figure_legend(fig, lang, M, ...
        'scatter', brine);
    add_mc_condition_header(fig, lang, cn_note, en_note, []);
    [ellipse_note, ellipse_interp, ellipse_font] = figure_label(choose_label(lang, ...
        '虚线椭圆：95% 高斯协方差范围', 'Dashed ellipses: 95% Gaussian covariance regions'), lang);
    annotation(fig, 'textbox', [.16 .008 .80 .043], 'String', ellipse_note, ...
        'Interpreter', ellipse_interp, 'FontName', ellipse_font, 'FontSize', 8.5, ...
        'EdgeColor', 'none', 'Margin', 0, 'HorizontalAlignment', 'center');
else
    Pe = draw_mc_lda_panel(fig, M, lang, Ip_c, R_c, Ip_h, R_h, stat);
    add_mc_figure_legend(fig, lang, M, 'lda', false);
    add_mc_condition_header(fig, lang, cn_note, en_note, Pe);
end
format_landscape_mc_panel(fig);
add_panel_id(fig, ['(' id ')']);
assert_no_panel_titles(fig);
end

function draw_mc_scatter_panel(fig, M, lang, Ip_c, R_c, Ip_h, R_h, show_brine)
xl_all = [Ip_c; Ip_h]; yl_all = [R_c; R_h];
span_ip = max(max(xl_all)-min(xl_all), 1e-6);
span_r  = max(max(yl_all)-min(yl_all), 1e-6);
xl = [min(xl_all)-.10*span_ip, max(xl_all)+.10*span_ip];
yl = [min(yl_all)-.10*span_r, max(yl_all)+.10*span_r];

% Compact but physically legible 6.40 x 3.60 in landscape layout.
% Leave a dedicated band above top KDE for conditions and legend.
main_x = .145; main_y = .185; main_w = .635; main_h = .530;
top_gap = .010; top_h = .064;
side_gap = .012; side_w = .063;
ax = axes('Parent', fig, 'Units', 'normalized', ...
    'Position', [main_x main_y main_w main_h], 'Color', 'white');
hold(ax, 'on');
if show_brine
    plot(ax, M.Ip_sw_line, M.VpVs_sw_line, '-', ...
        'Color', .76*M.col_sw+.24*[1 1 1], 'LineWidth', 2.5, ...
        'HandleVisibility', 'off');
end
nshow = min(1900, numel(Ip_c));
scatter(ax, Ip_c(1:nshow), R_c(1:nshow), 12, ...
    .72*M.col_ch4+.28*[1 1 1], 'filled', ...
    'MarkerEdgeColor', 'none', 'HandleVisibility', 'off');
scatter(ax, Ip_h(1:nshow), R_h(1:nshow), 12, ...
    .72*M.col_h2+.28*[1 1 1], 'filled', ...
    'MarkerEdgeColor', 'none', 'HandleVisibility', 'off');
add_covariance_ellipse(ax, [Ip_c R_c], M.col_ch4, .95);
add_covariance_ellipse(ax, [Ip_h R_h], M.col_h2, .95);
xlim(ax, xl); ylim(ax, yl);
style_axes(ax, 18.3);
apply_axis_labels(ax, lang, ...
    '纵波阻抗 Iₚ (MRayl)', 'P-wave impedance Iₚ (MRayl)', ...
    '纵横波速比 Vₚ/Vₛ', 'Vₚ/Vₛ');
if span_ip < .25
    step = nice_tick_step(diff(xl)/5);
    t0 = ceil((xl(1)-1e-10)/step)*step;
    t1 = floor((xl(2)+1e-10)/step)*step;
    xticks(ax, t0:step:t1);
    xtickformat(ax, '%.2f');
end
if span_r < .04, ytickformat(ax, '%.3f'); end
grid(ax, 'on');

% Separate marginal density axes have no ticks or labels, and no title.
ax_top = axes('Parent', fig, 'Units', 'normalized', ...
    'Position', [main_x main_y+main_h+top_gap main_w top_h], ...
    'Color', 'white');
hold(ax_top, 'on');
xi = linspace(xl(1), xl(2), 220);
draw_reference_marginal(ax_top,Ip_c,Ip_h,xi,M.col_ch4,M.col_h2,'top');

ax_side = axes('Parent', fig, 'Units', 'normalized', ...
    'Position', [main_x+main_w+side_gap main_y side_w main_h], ...
    'Color', 'white');
hold(ax_side, 'on');
yi = linspace(yl(1), yl(2), 220);
draw_reference_marginal(ax_side,R_c,R_h,yi,M.col_ch4,M.col_h2,'side');
end

function Pe = draw_mc_lda_panel(fig,M,lang,Ip_c,R_c,Ip_h,R_h,stat) %#ok<INUSD>
% Display the same stable LDA direction and KDE overlap used by the maps.
% Point masses have probability mass markers, not an invented noise width.
Pe=stat.Pe_kde;
ax=axes('Parent',fig,'Units','normalized', ...
    'Position',[.145 .185 .735 .594],'Color','white');
hold(ax,'on'); style_axes(ax,18.3);
if isempty(stat.z_c) || isempty(stat.z_h)
    [note,it,fnt]=figure_label(choose_label(lang, ...
        'LDA 均值方向未定义；请查看二维高斯结果', ...
        'LDA mean direction is undefined; use the 2D Gaussian result'),lang);
    text(ax,.5,.5,note,'Units','normalized','Interpreter',it,'FontName',fnt, ...
        'HorizontalAlignment','center','FontSize',11);
    set(ax,'XTick',[],'YTick',[]);
    return;
end
[z,pdf_c,pdf_h,atoms]=reference_projection_densities(stat);
colors={M.col_ch4,M.col_h2}; loc=[stat.z_c(1),stat.z_h(1)];
if all(atoms)
    for k=1:2
        stem(ax,loc(k),1,'Color',colors{k},'MarkerFaceColor',colors{k}, ...
            'LineWidth',2,'HandleVisibility','off');
    end
    if loc(1)==loc(2)
        stem(ax,loc(1),1,'Color',[.80 .18 .25],'Marker','s', ...
            'LineWidth',1.3,'HandleVisibility','off');
    end
    span=max(diff([min(loc) max(loc)]),1);
    xlim(ax,[min(loc)-.2*span max(loc)+.2*span]); ylim(ax,[0 1.1]);
    apply_axis_labels(ax,lang, ...
        '标准化 LDA 投影得分 η','Standardized LDA projection score η', ...
        '概率质量','Probability mass');
else
    if ~atoms(1)
        fill(ax,z,pdf_c,.20*M.col_ch4+.80*[1 1 1], ...
            'EdgeColor',M.col_ch4,'LineWidth',2.25,'HandleVisibility','off');
    end
    if ~atoms(2)
        fill(ax,z,pdf_h,.20*M.col_h2+.80*[1 1 1], ...
            'EdgeColor',M.col_h2,'LineWidth',2.25,'HandleVisibility','off');
    end
    fill(ax,z,min(pdf_c,pdf_h),[.93 .75 .77], ...
        'EdgeColor',[.80 .18 .25],'LineWidth',1.3,'HandleVisibility','off');
    for k=find(atoms)
        xline(ax,loc(k),'-','Color',colors{k},'LineWidth',2,'HandleVisibility','off');
    end
    peak=max([pdf_c(:);pdf_h(:)]);
    ylim(ax,[0 1.07*max(peak,realmin)]); xlim(ax,[z(1) z(end)]);
    apply_axis_labels(ax,lang, ...
        '标准化 LDA 投影得分 η','Standardized LDA projection score η', ...
        '概率密度','Probability density');
    if any(atoms)
        [note,it,fnt]=figure_label(choose_label(lang, ...
            '竖线表示点质量','Vertical line: point mass'),lang);
        text(ax,.98,.95,note,'Units','normalized','Interpreter',it,'FontName',fnt, ...
            'FontSize',9,'HorizontalAlignment','right','BackgroundColor','white');
    end
end
xline(ax,0,':','LineWidth',1.1,'Color',[.32 .36 .41],'HandleVisibility','off');
grid(ax,'on');
end

function step = nice_tick_step(raw_step)
% 1/2/5/10-based intervals; keeps (c) ticks readable and familiar.
exponent = floor(log10(max(raw_step, realmin)));
scale = 10^exponent;
r = raw_step/scale;
if r <= 1, b = 1;
elseif r <= 2, b = 2;
elseif r <= 2.6, b = 2;
elseif r <= 5, b = 5;
else, b = 10;
end
step = b*scale;
end

function add_mc_figure_legend(fig, lang, M, kind, show_brine)
% Independent 6.40 x 3.60-in panel: one compact, readable legend per plot.
% Legend occupies the top strip; no data axes or plot titles are covered.
y = .952;
if strcmp(kind, 'scatter')
    if show_brine
        annotation(fig, 'line', [.095 .133], [y y], ...
            'Color', M.col_sw, 'LineWidth', 1.9);
        mc_legend_text(fig, lang, .143, y, ...
            '盐水参考', 'Brine reference', .195);
        ch4_x = .385;
        h2_x = .675;
    else
        ch4_x = .245;
        h2_x = .565;
    end
    annotation(fig, 'ellipse', [ch4_x y-.009 .014 .018], ...
        'FaceColor', M.col_ch4, 'EdgeColor', M.col_ch4);
    mc_legend_text(fig, lang, ch4_x+.022, y, ...
        '甲烷 (CH₄)', 'Methane (CH₄)', .225);
    annotation(fig, 'ellipse', [h2_x y-.009 .014 .018], ...
        'FaceColor', M.col_h2, 'EdgeColor', M.col_h2);
    mc_legend_text(fig, lang, h2_x+.022, y, ...
        '氢气 (H₂)', 'Hydrogen (H₂)', .215);
else
    annotation(fig, 'ellipse', [.105 y-.009 .014 .018], ...
        'FaceColor', M.col_ch4, 'EdgeColor', M.col_ch4);
    mc_legend_text(fig, lang, .128, y, ...
        '甲烷 (CH₄)', 'Methane (CH₄)', .215);
    annotation(fig, 'ellipse', [.375 y-.009 .014 .018], ...
        'FaceColor', M.col_h2, 'EdgeColor', M.col_h2);
    mc_legend_text(fig, lang, .398, y, ...
        '氢气 (H₂)', 'Hydrogen (H₂)', .205);
    annotation(fig, 'rectangle', [.650 y-.008 .020 .016], ...
        'FaceColor', [.93 .75 .77], 'EdgeColor', [.80 .18 .25], ...
        'LineWidth', 0.9);
    mc_legend_text(fig, lang, .680, y, ...
        'LDA 重叠部分', 'LDA overlap', .275);
end
annotation(fig, 'line', [.075 .960], [.908 .908], ...
    'Color', [.87 .89 .91], 'LineWidth', .60);
end

function mc_legend_text(fig, lang, x, y, cn, en, w)
[legend_text, legend_intrp, legend_font] = figure_label( ...
    choose_label(lang, cn, en), lang);
annotation(fig, 'textbox', [x, y-.028, w, .056], ...
    'String', legend_text, 'Interpreter', legend_intrp, ...
    'FontName', legend_font, 'FontSize', 10.3, ...
    'FontWeight', 'bold', 'EdgeColor', 'none', ...
    'Color', [.19 .22 .26], ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'FitBoxToText', 'off', 'Margin', 0);
end

function add_mc_condition_header(fig, lang, cn_note, en_note, Pe)
% Band below legend and above KDE/main axes; (a)-(d) sits at far left.
[condition, intrp, fnt] = figure_label( ...
    choose_label(lang, cn_note, en_note), lang);
if isempty(Pe)
    note_w = .790;
else
    note_w = .600;  % reserve right edge for the LDA overlap result
end
annotation(fig, 'textbox', [.160 .815 note_w .068], ...
    'String', condition, 'Interpreter', intrp, 'FontName', fnt, ...
    'FontSize', 10.5, 'FontWeight', 'bold', ...
    'EdgeColor', 'none', 'Color', [.17 .19 .22], ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'FitBoxToText', 'off', 'Margin', 0);
if ~isempty(Pe)
    if isfinite(Pe), raw_pe = sprintf('P_e (KDE) = %.2f%%',Pe*100);
    else, raw_pe = 'P_e (KDE) = n/a'; end
    [pe_text, pe_interp] = scientific_label(raw_pe);
    annotation(fig, 'textbox', [.739 .815 .223 .068], ...
        'String', pe_text, 'Interpreter', pe_interp, ...
        'FontName', 'Times New Roman', 'FontSize', 10.5, ...
        'FontWeight', 'bold', 'EdgeColor', 'none', ...
        'Color', [.17 .19 .22], 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'middle', ...
        'FitBoxToText', 'off', 'Margin', 0);
end
end

function format_landscape_mc_panel(fig)
% Same effective PDF font sizes as Figure 1(a,b,d,e,f,g):
% 10.3-pt ticks and 11.5-pt visible x/y labels (TeX has export correction).
page = getappdata(fig, 'PhysicalExportSizeInches');
assert(is_landscape_panel_size(page), 'RPT:InvalidMCLandscape', ...
    'Every Figure 2 subplot must be 6.40 x 3.60 inches.');
axs = findall(fig, 'Type', 'Axes');
for k = 1:numel(axs)
    ax = axs(k);
    if strcmp(get(ax, 'Visible'), 'off')
        continue; % marginal KDE axes intentionally have no tick labels
    end
    set(ax, 'FontSize', 10.3, 'LineWidth', 1.15, ...
        'TickLength', [.015 .018]);
    set(get(ax, 'XLabel'), 'FontUnits', 'points', ...
        'FontSize', axis_label_font_size(get(get(ax, 'XLabel'), 'Interpreter')), ...
        'FontWeight', 'bold');
    set(get(ax, 'YLabel'), 'FontUnits', 'points', ...
        'FontSize', axis_label_font_size(get(get(ax, 'YLabel'), 'Interpreter')), ...
        'FontWeight', 'bold');
end
end

function lock_final_panel_typography(fig)
% Last-pass export safeguard.  FontSize alone is not enough: TeX and
% non-TeX labels with numerically equal font sizes differ in actual PDF.
axs = findall(fig, 'Type', 'Axes');
for k = 1:numel(axs)
    ax = axs(k);
    if strcmp(get(ax, 'Visible'), 'off'), continue; end
    set(ax, 'FontUnits', 'points', 'FontSize', 10.3);
    hx = get(ax, 'XLabel');
    hy = get(ax, 'YLabel');
    set(hx, 'FontUnits', 'points', ...
        'FontSize', axis_label_font_size(get(hx, 'Interpreter')), ...
        'FontWeight', 'bold');
    set(hy, 'FontUnits', 'points', ...
        'FontSize', axis_label_font_size(get(hy, 'Interpreter')), ...
        'FontWeight', 'bold');
end
cb = findall(fig, 'Type', 'ColorBar');
for k = 1:numel(cb)
    set(cb(k), 'FontSize', 10.3); % ColorBar.FontSize is already in points.
    lab = get(cb(k), 'Label');
    set(lab, 'FontUnits', 'points', ...
        'FontSize', colorbar_label_font_size(get(lab, 'Interpreter')), ...
        'FontWeight', 'bold');
end
% Fix physical axes coordinates AFTER font/colorbar layout changes.
% All 6 portrait panels must have the exact same 1.392 x 2.520 in plot box.
page_in = getappdata(fig, 'PhysicalExportSizeInches');
if page_in(1) < page_in(2)
    main_ax = axs(arrayfun(@(h) strcmp(get(h,'Visible'),'on'), axs));
    assert(numel(main_ax)==1, 'RPT:PortraitAxisCount', ...
        'Every portrait panel must have exactly one main axis.');
    set(main_ax, 'Units', 'normalized', ...
        'Position', [.235 .150 .480 .700]);
    if numel(cb)==1
        set(cb, 'Units', 'normalized', ...
            'Position', [.766 .150 .038 .700]);
        % MATLAB colorbar auto-layout may move axes after positioning the
        % bar; explicitly restore the fixed data rectangle one last time.
        set(main_ax, 'Position', [.235 .150 .480 .700]);
    end
end
end

function assert_consistent_main_axis_fonts(fig)
% Prevent future edits from introducing Figure 1 / Figure 2 font drift.
% Marginal KDE axes in Figure 2(a,c) are invisible and have no text.
axs = findall(fig, 'Type', 'Axes');
seen_main = 0;
for k = 1:numel(axs)
    ax = axs(k);
    if strcmp(get(ax, 'Visible'), 'off'), continue; end
    seen_main = seen_main + 1;
    assert(abs(get(ax, 'FontSize') - 10.3) < 0.01, ...
        'RPT:InconsistentTickFont', 'Main axes ticks must be 10.3 pt.');
    hx = get(ax, 'XLabel'); hy = get(ax, 'YLabel');
    fs_x = axis_label_font_size(get(hx, 'Interpreter'));
    fs_y = axis_label_font_size(get(hy, 'Interpreter'));
    assert(strcmpi(get(hx, 'FontUnits'), 'points') && ...
           strcmpi(get(hy, 'FontUnits'), 'points') && ...
           abs(get(hx, 'FontSize') - fs_x) < 0.01 && ...
           abs(get(hy, 'FontSize') - fs_y) < 0.01, ...
        'RPT:InconsistentAxisFont', ...
        'Main X/Y axis labels must have the same effective PDF size as Fig 1(a).');
end
assert(seen_main == 1, 'RPT:MultipleMainAxes', ...
    'Each standalone panel must contain exactly one visible main axes.');
% Enforce identical main-plot physical geometry for Figure 1's six
% portrait panels. Portrait output pages already have one fixed size.
page_in = getappdata(fig, 'PhysicalExportSizeInches');
if page_in(1) < page_in(2)
    main_ax = axs(arrayfun(@(h) strcmp(get(h,'Visible'),'on'), axs));
    old_units = get(main_ax, 'Units');
    set(main_ax, 'Units', 'normalized');
    pos = get(main_ax, 'Position');
    set(main_ax, 'Units', old_units);
    expected = [.235 .150 .480 .700];
    assert(all(abs(pos-expected) < 0.025), ...
        'RPT:PortraitGeometryDrift', ...
        'All 6 Figure 1 portrait main axes must share one physical rectangle.');
end
end

function draw_vector_heatmap(ax, X, Y, Z)
% Fixed 48-level filled contours are vector patches (not bitmap images).
% Numerical results and color limits remain unchanged.
lo = min(Z(:)); hi = max(Z(:));
if hi > lo
    levels = linspace(lo-1e-9*max(1,abs(lo)), ...
        hi+1e-9*max(1,abs(hi)), 49);
    contourf(ax, X, Y, Z, levels, 'LineStyle', 'none');
else
    surf(ax, X, Y, zeros(size(Z)), Z, 'FaceColor', 'flat', ...
        'EdgeColor', 'none');
    view(ax, 2);
end
set(ax, 'YDir', 'normal');
end

function add_small_condition_note(ax, lines, corner)
if strcmp(corner, 'southeast')
    xy = [.978 .045]; va = 'bottom';
else
    xy = [.978 .955]; va = 'top';
end
[display_lines, intrp] = scientific_label(lines);
text(ax, xy(1), xy(2), display_lines, ...
    'Units', 'normalized', 'Interpreter', intrp, ...
    'FontName', 'Times New Roman', 'FontSize', 14.8, ...
    'FontWeight', 'bold', 'Color', [.16 .19 .23], ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', va, ...
    'BackgroundColor', 'w', 'Margin', 1.6, 'Clipping', 'on');
end

function draw_rpt_key_inside(ax, lang, V, col_sw, col_ch4, col_h2)
% Figure 1(c): five small keyed strokes/markers occupy a permanently
% reserved EMPTY STRIP at the BOTTOM of the RPT axes. This approach is
% independent of tiledlayout's inter-row spacing, figure/label font
% metrics and print/pdf renderer legend-box inflation.
%
% The lowest RPT grid point is forced above 21.5% of the plotted height.
% The key's upper text edge is ~15%; it cannot obscure the actual mesh.
mesh_min = inf;
for j = 1:numel(V.Sat_C)
    mesh_min = min(mesh_min, min(V.Sat_C{j}(:,2)));
end
for j = 2:numel(V.Sat_H)
    mesh_min = min(mesh_min, min(V.Sat_H{j}(:,2)));
end
for j = 1:numel(V.Phi_C)
    mesh_min = min(mesh_min, min(V.Phi_C{j}(:,2)));
    mesh_min = min(mesh_min, min(V.Phi_H{j}(:,2)));
end
yl = ylim(ax);
target = .255;
if (mesh_min - yl(1)) / diff(yl) < target
    new_ymin = (mesh_min - target*yl(2))/(1-target);
    ylim(ax, [new_ymin yl(2)]);
end
% Only the three curve-family keys are shown. Reference markers remain
% in the RPT data, but their explanatory text entries are intentionally omitted.
draw_rpt_key_entry(ax, lang, .047, .142, 'line', ...
    col_sw,  '盐水',       'Brine');
draw_rpt_key_entry(ax, lang, .358, .142, 'dash', ...
    col_ch4, '甲烷网格',   'CH₄ grid');
draw_rpt_key_entry(ax, lang, .667, .142, 'line', ...
    col_h2,  '氢气网格',   'H₂ grid');
end

function draw_rpt_key_entry(ax, lang, u, v, kind, color, cn, en)
% All symbols and text use the axes coordinate range, so none can extend
% across the xlabel or spill into the next row.
xl = xlim(ax); yl = ylim(ax);
dx = diff(xl); dy = diff(yl);
x0 = xl(1) + u*dx; y0 = yl(1) + v*dy;
text_x = x0 + .055*dx;
if strcmp(kind, 'line') || strcmp(kind, 'dash')
    ls = '-';
    if strcmp(kind, 'dash'), ls = '--'; end
    plot(ax, [x0 x0+.043*dx], [y0 y0], ls, ...
        'Color', color, 'LineWidth', 2.25, ...
        'Clipping', 'on', 'HandleVisibility', 'off');
else
    mk = 'o';
    if strcmp(kind, 'diamond'), mk = 'd'; end
    plot(ax, x0+.019*dx, y0, mk, 'Color', color, ...
        'MarkerFaceColor', color, 'MarkerEdgeColor', [.20 .23 .27], ...
        'MarkerSize', 6.0, 'LineWidth', .9, ...
        'Clipping', 'on', 'HandleVisibility', 'off');
end
[item_text, intrp, item_font] = figure_label( ...
    choose_label(lang, cn, en), lang);
text(ax, text_x, y0, item_text, ...
    'Interpreter', intrp, 'FontName', item_font, ...
    'FontSize', 10.4, 'FontWeight', 'bold', ...
    'Color', [.19 .22 .26], 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'middle', 'Clipping', 'on', ...
    'HandleVisibility', 'off');
end

function style_axes(ax, sz)
set(ax, 'FontName', 'Times New Roman', 'FontSize', sz, ...
    'FontWeight', 'bold', 'LineWidth', 1.95, ...
    'TickDir', 'out', 'TickLength', [.020 .022], ...
    'Box', 'on', 'Layer', 'top', ...
    'XColor', [.18 .22 .26], 'YColor', [.18 .22 .26], ...
    'GridLineStyle', '-', 'GridColor', [.35 .39 .44], ...
    'GridAlpha', .15, 'MinorGridAlpha', .10);
end

function apply_axis_labels(ax, lang, cn_x, en_x, cn_y, en_y)
% figure_label() ensures Chinese labels are truly mixed-font: Chinese
% glyphs in a detected CJK font; Latin, Greek, subscripts and units in Times.
[x_text, x_interpreter, x_font] = figure_label( ...
    choose_label(lang, cn_x, en_x), lang);
[y_text, y_interpreter, y_font] = figure_label( ...
    choose_label(lang, cn_y, en_y), lang);
xlabel(ax, x_text, ...
    'Interpreter', x_interpreter, 'FontName', x_font, ...
    'FontUnits', 'points', 'FontSize', axis_label_font_size(x_interpreter), ...
    'FontWeight', 'bold', 'Color', [.11 .15 .21]);
ylabel(ax, y_text, ...
    'Interpreter', y_interpreter, 'FontName', y_font, ...
    'FontUnits', 'points', 'FontSize', axis_label_font_size(y_interpreter), ...
    'FontWeight', 'bold', 'Color', [.11 .15 .21]);
end

function fs = axis_label_font_size(interpreter)
% Target physical, visible font size for every main X/Y axis label: 11.5 pt.
% The supplied Figure 1(a) EN PDF shows "Porosity" (Interpreter='none')
% printed at ~11.4 pt, while "Gas saturation S_{g}" (Interpreter='tex')
% prints at ~9.0 pt despite both being set to FontSize=11.5 in MATLAB.
% The MATLAB TeX/PDF renderer reduces its effective PDF text size.  Apply
% an empirical TeX-only correction; do NOT change the scientific symbols
% or convert subscripts to Unicode, which causes font substitution.
if nargin < 1 || isempty(interpreter), interpreter = 'none'; end
fs = 11.5;
if strcmpi(char(interpreter), 'tex')
    fs = fs * (11.4 / 9.0);  % ~14.57 MATLAB pt -> ~11.5 pt in print
end
end

function fs = colorbar_label_font_size(interpreter)
% Physical colorbar label size also matches the Fig 1(a) X-axis label.
% In the supplied PDF, TeX colorbar label at MATLAB 10.0 pt prints ~8.4 pt.
if nargin < 1 || isempty(interpreter), interpreter = 'none'; end
fs = 11.5;
if strcmpi(char(interpreter), 'tex')
    fs = fs * (10.0 / 8.4);  % TeX/colorbar PDF-specific correction
end
end

function style_colorbar(ax, lang, cn_text, en_text)
cb = colorbar(ax);
set(cb, 'FontName', 'Times New Roman', 'FontSize', 10.3, ...
    'FontWeight', 'bold', 'LineWidth', 1.45, 'TickDirection', 'out');
[cb_label, cb_interp, cb_font] = figure_label( ...
    choose_label(lang, cn_text, en_text), lang);
cb.Label.String = cb_label;
cb.Label.Interpreter = cb_interp;
cb.Label.FontName = cb_font;
cb.Label.FontUnits = 'points';
cb.Label.FontSize = colorbar_label_font_size(cb_interp);
cb.Label.FontWeight = 'bold';
cb.Label.Color = [.11 .15 .21];
end

function add_panel_id(fig, panel_id)
% Keep (a)-(g)/(a)-(d), but NEVER add subplot names.
page_in = getappdata(fig, 'PhysicalExportSizeInches');
if ~isempty(page_in) && page_in(1) < page_in(2)
    % Small 3.6-inch portrait panel: letter in its own narrow top band.
    pos = [.075 .931 .13 .046];
    fontsize = 14.2;
else
    % 6.40 x 3.60-in landscape: ID in its own left-hand header gutter.
    pos = [.048 .815 .080 .068];
    fontsize = 14.2;
end
annotation(fig, 'textbox', pos, ...
    'String', panel_id, 'Interpreter', 'none', ...
    'FontName', 'Times New Roman', 'FontSize', fontsize, ...
    'FontWeight', 'bold', 'Color', [.10 .13 .17], ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'middle', 'FitBoxToText', 'off', 'Margin', 0);
end

function assert_no_panel_titles(fig)
% Fail fast if a future revision accidentally calls title() or sgtitle().
axs = findall(fig, 'Type', 'Axes');
for k = 1:numel(axs)
    t = get(axs(k), 'Title');
    if ~isempty(strtrim(char(string(get(t, 'String')))))
        error('RPT:UnexpectedPanelTitle', ...
            'Panel contains a subplot name/title; remove title() before export.');
    end
end
if ~isempty(findall(fig, 'Type', 'Legend'))
    error('RPT:UnexpectedLegend', ...
        'Panel legends must use fixed annotations, not auto-positioned legends.');
end
end

function add_condition_note(ax, content, corner)
if strcmp(corner, 'northeast')
    pos = [.98 .975]; horiz = 'right';
else
    pos = [.025 .975]; horiz = 'left';
end
[note_text, intrp] = scientific_label(content);
text(ax, pos(1), pos(2), note_text, ...
    'Units', 'normalized', 'Interpreter', intrp, ...
    'FontName', 'Times New Roman', 'FontSize', 16.0, ...
    'FontWeight', 'bold', 'Color', [.17 .21 .25], ...
    'HorizontalAlignment', horiz, 'VerticalAlignment', 'top', ...
    'BackgroundColor', [1 1 1], 'Margin', 4, 'Clipping', 'on');
end

function [out, interpreter, font] = figure_label(raw, lang)
% FIGURE_LABEL is the ONLY entry point for localized figure labels.
% MATLAB Text objects have one FontName, so assigning SimSun to the entire
% Chinese label forces Latin letters and numerical units into the CJK face.
% This function explicitly switches fonts within a mixed Chinese/Latin
% string. For labels without Chinese text there is NO font switch and all
% scientific glyphs are rendered with Times New Roman.
[out, interpreter] = scientific_label(raw);
font = 'Times New Roman';
if ~strcmp(lang, 'CN'), return; end
cn_font = chinese_figure_font();  % resolved once and cached for this run
if iscell(out)
    % Not used by our present captions but supported for future additions.
    tmp = cell(size(out));
    for k = 1:numel(out)
        [tmp{k}, intrp_k] = figure_label(out{k}, lang);
        if strcmp(intrp_k, 'tex'), interpreter = 'tex'; end
    end
    out = tmp;
    return;
end

% Cover common Han characters, CJK punctuation and full-width symbols.
c = double(out);
is_cn = (c >= 13312 & c <= 40959) | (c >= 63744 & c <= 64255) | ...
        (c >= 12288 & c <= 12351) | (c >= 65280 & c <= 65519);
if ~any(is_cn), return; end

% Pure Chinese + punctuation uses the selected CJK font.
if ~any(~is_cn & ~isspace(out))
    font = cn_font;
    interpreter = 'none';
    return;
end

% A genuine mixed label: keep Chinese runs in CJK and Latin/numeric
% runs Times. TeX subscripts are grouped for the CHINESE mixed-font case,
% which prevents a Tahoma glyph fallback for Unicode subscripts p/s/2/4.
% English-only labels stay in searchable Unicode and Times New Roman.
out = strrep(out, 'ₚ', '_{p}');
out = strrep(out, 'ₛ', '_{s}');
out = strrep(out, 'ₑ', '_{e}');
out = strrep(out, '₄', '_{4}');
out = strrep(out, '₂', '_{2}');
% Font commands and TeX subscript syntax contain no CJK glyphs, so we must
% update the character mask AFTER this expansion.
c = double(out);
is_cn = (c >= 13312 & c <= 40959) | (c >= 63744 & c <= 64255) | ...
        (c >= 12288 & c <= 12351) | (c >= 65280 & c <= 65519);
buf = '\bf'; % TeX uses its own bold switch; FontWeight alone may not apply.
last_cn = false;
for k = 1:numel(out)
    now_cn = is_cn(k);
    if k == 1 || (now_cn ~= last_cn)
        if now_cn
            buf = [buf '\fontname{' cn_font '}']; %#ok<AGROW>
        else
            buf = [buf '\fontname{Times New Roman}']; %#ok<AGROW>
        end
    end
    buf = [buf out(k)]; %#ok<AGROW>
    last_cn = now_cn;
end
out = buf;
interpreter = 'tex';
end

function s = choose_label(lang, cn_text, en_text)
if strcmp(lang, 'CN')
    s = cn_text;
else
    s = en_text;
end
end

function [out, interpreter] = scientific_label(raw)
% SCIENTIFIC_LABEL: normalize scientific notation for vector PDF/SVG export.
% All parameter subscripts are unified as TeX grouped subscripts so that
% I_p, V_p/V_s, S_g, P_e, CH_4 and H_2 keep correct baseline alignment.
if iscell(raw)
    out = cell(size(raw));
    interpreter = 'none';
    for j = 1:numel(raw)
        [out{j}, it] = scientific_label(raw{j});
        if strcmp(it, 'tex'), interpreter = 'tex'; end
    end
    return;
end
if isstring(raw)
    assert(isscalar(raw), 'RPT:NonScalarLabel', ...
        'A scientific text string must be scalar.');
    raw = char(raw);
end
assert(ischar(raw), 'RPT:InvalidLabelType', ...
    'Expected a char/string or a cell array of char labels.');

out = raw;
out = strrep(out, '\phi', 'φ');
out = strrep(out, '\eta', 'η');
out = strrep(out, '^{\circ}', '°');
out = strrep(out, '\circ', '°');

% Accept both ASCII-style and Unicode-subscript spellings, then unify all
% of them to grouped TeX subscripts.
subs = {
    'I_p',   'I_{p}';
    'Iₚ',    'I_{p}';
    'V_p',   'V_{p}';
    'Vₚ',    'V_{p}';
    'V_s',   'V_{s}';
    'Vₛ',    'V_{s}';
    'P_e',   'P_{e}';
    'Pₑ',    'P_{e}';
    'S_g',   'S_{g}';
    'CH_4',  'CH_{4}';
    'CH₄',   'CH_{4}';
    'H_2',   'H_{2}';
    'H₂',    'H_{2}';
};
for k = 1:size(subs,1)
    out = strrep(out, subs{k,1}, subs{k,2});
end

if contains(out, '_{')
    interpreter = 'tex';
else
    interpreter = 'none';
end

if ~isempty(strfind(out, '\phi')) || ...
        ~isempty(strfind(out, '\eta')) || ...
        ~isempty(strfind(out, '^{\circ}'))
    error('RPT:TypographyError', 'An unnormalized scientific label remains.');
end
end

function check_figure_fonts()
% SimSun is OPTIONAL. MATLAB installs on macOS/Linux often do not include it.
% Times New Roman is kept strict as it is the requested scientific font.
try
    available = listfonts;
catch fontErr
    warning('RPT:FontListUnavailable', ...
        'Could not query installed fonts: %s', fontErr.message);
    available = {};
end
if ischar(available), available = cellstr(available); end
if isstring(available), available = cellstr(available); end
if ~isempty(available) && ~any(strcmpi(available, 'Times New Roman'))
    error('RPT:TimesFontUnavailable', ...
        'Times New Roman is required for scientific labels/units. Install it, then rerun.');
end
cn_font = chinese_figure_font();
fprintf('Figure typography: scientific = Times New Roman; Chinese = %s.\n', cn_font);
end

function font = chinese_figure_font()
% Return the REAL font name reported by listfonts, not an assumed alias.
% Prefer serif faces similar to SimSun, then readable sans-serif CJK faces.
% No hard failure when SimSun itself is absent.
persistent cached_cn_font
if ~isempty(cached_cn_font)
    font = cached_cn_font;
    return;
end

try
    installed = listfonts;
    if ischar(installed), installed = cellstr(installed); end
    if isstring(installed), installed = cellstr(installed); end
catch fontErr
    warning('RPT:CJKFontCheckSkipped', ...
        'Could not inspect Chinese fonts (%s). Trying MATLAB default.', fontErr.message);
    installed = {};
end

preferred = { ...
    'SimSun', 'NSimSun', 'Songti SC', 'Songti TC', 'STSong', ...
    'Noto Serif CJK SC', 'Noto Serif CJK TC', 'Noto Serif SC', ...
    'Source Han Serif SC', 'Source Han Serif CN', ...
    'Microsoft YaHei', 'Microsoft YaHei UI', 'PingFang SC', ...
    'Noto Sans CJK SC', 'Noto Sans SC', ...
    'Source Han Sans SC', 'Source Han Sans CN', ...
    'SimHei', 'Heiti SC', 'STHeiti', 'DengXian', ...
    'WenQuanYi Micro Hei', 'WenQuanYi Zen Hei', ...
    'Arial Unicode MS', 'KaiTi', 'FangSong'};
font = '';
for k = 1:numel(preferred)
    found = find(strcmpi(installed, preferred{k}), 1, 'first');
    if ~isempty(found)
        font = installed{found};
        break;
    end
end

% Some installations use versioned/regional font family names.
if isempty(font)
    terms = {'CJK', 'Source Han', 'YaHei', 'PingFang', 'Songti', ...
             'SimHei', 'SimSun', 'Noto Sans SC', 'Noto Serif SC', ...
             'WenQuanYi', 'Heiti', 'DengXian', 'KaiTi', 'FangSong'};
    for k = 1:numel(terms)
        found = find(contains(installed, terms{k}, 'IgnoreCase', true), 1, 'first');
        if ~isempty(found)
            font = installed{found};
            break;
        end
    end
end

if isempty(font)
    % Do not abort an otherwise valid numerical script. MATLAB's font
    % substitution may still render Han glyphs with its platform fallback.
    font = get(groot, 'DefaultTextFontName');
    if isstring(font), font = char(font); end
    if isempty(font), font = 'Helvetica'; end
    warning('RPT:NoRecognizedCJKFont', ...
        ['No recognized Chinese font found. Using MATLAB default (%s). ' ...
         'Chinese characters may be substituted or missing in exported PDFs. ' ...
         'Install a CJK font such as Noto Serif CJK SC if needed.'], font);
elseif ~strcmpi(font, 'SimSun')
    fprintf('SimSun not found; using installed Chinese font: %s.\n', font);
end
cached_cn_font = font;
end

function validate_figure_typography(fig)
% Fail before export if any legacy scientific-token spellings remain in
% visible vector text objects. The PDF thus comes from normalized labels.
labels = findall(fig, 'Type', 'Text');
for j = 1:numel(labels)
    str = get(labels(j), 'String');
    if iscell(str)
        str = sprintf('%s ', str{:});
    elseif isstring(str)
        str = char(join(str, ' '));
    end
    if ~ischar(str), continue; end
    forbidden = {'\phi','\eta','^{\circ}', 'I_p','V_p','V_s', ...
                 'P_e','CH_4','H_2', 'S_g'};
    for k = 1:numel(forbidden)
        if ~isempty(strfind(str, forbidden{k})) %#ok<STREMP>
            error('RPT:LegacyTypography', ...
                'The figure contains legacy text markup: %s', forbidden{k});
        end
    end
end
end

function [C_joint, C_ip, C_ratio] = calc_ideal_contrast(Ip_h, R_h, Ip_c, R_c)
% 理想模型的 H2-CH4 相对弹性差异 [百分比]。
% 取两种气体预测值的平均作为参考值，避免偏向任何一种气体。
% 两个无量纲的百分比属性等权 RMS 组合，与误差统计无关。
Ip_ref = 0.5*(Ip_h + Ip_c);
R_ref  = 0.5*(R_h + R_c);
C_ip   = 100*abs(Ip_h - Ip_c) ./ max(abs(Ip_ref), realmin);
C_ratio = 100*abs(R_h - R_c) ./ max(abs(R_ref), realmin);
C_joint = sqrt(0.5*(C_ip.^2 + C_ratio.^2));
end

function pdf_vals = calc_1d_kde(samples, x_grid)
samples = samples(:);
n = length(samples);
bw = 1.06 * std(samples) * (n^(-0.2));
if bw <= 0, bw = 1e-4; end
diff_mat = (x_grid(:)' - samples) / bw;
pdf_vals = mean(exp(-0.5 * diff_mat.^2), 1) / (bw * sqrt(2*pi));
end

function lut = build_coolprop_lut(T_vec, P_vec)
[TT, PP] = meshgrid(T_vec, P_vec);
[rho_c, K_c] = get_gas_props_vec('Methane',  TT, PP);
[rho_h, K_h] = get_gas_props_vec('Hydrogen', TT, PP);
lut.T_vec = T_vec; lut.P_vec = P_vec;
lut.rho_c = rho_c; lut.K_c   = K_c;
lut.rho_h = rho_h; lut.K_h   = K_h;
end

function [Ip, VpVs] = calc_rpt_lut(Kd, Gd, phi, Sg, T_val, P_val, gasType, salinity, K0, rho0, lut)
T_grid = T_val + zeros(size(phi));
P_grid = P_val + zeros(size(phi));
if any(~isfinite(T_grid(:))) || any(~isfinite(P_grid(:))) || ...
        any(T_grid(:) < lut.T_vec(1) | T_grid(:) > lut.T_vec(end)) || ...
        any(P_grid(:) < lut.P_vec(1) | P_grid(:) > lut.P_vec(end))
    error('RPT:LUTOutOfRange', ...
        'Requested T/P lies outside the CoolProp LUT; extend the LUT range.');
end
if strcmpi(gasType, 'Methane')
    rho_g = interp2(lut.T_vec, lut.P_vec, lut.rho_c, T_grid, P_grid, 'spline');
    K_g   = interp2(lut.T_vec, lut.P_vec, lut.K_c,   T_grid, P_grid, 'spline');
else
    rho_g = interp2(lut.T_vec, lut.P_vec, lut.rho_h, T_grid, P_grid, 'spline');
    K_g   = interp2(lut.T_vec, lut.P_vec, lut.K_h,   T_grid, P_grid, 'spline');
end
[rho_w, K_w] = get_brine_props(T_grid, P_grid, salinity);
Sw = 1 - Sg;
Kbrine = gassmann_vec(Kd, K0, phi, K_w);
Kgas   = gassmann_vec(Kd, K0, phi, K_g);
c      = (4/3) .* Gd;
Ksat   = 1 ./ (Sw ./ (Kbrine + c) + Sg ./ (Kgas + c)) - c;
rho_fl = Sg .* rho_g + Sw .* rho_w;
rho_sat = (1 - phi) .* rho0 + phi .* rho_fl;
Gsat    = Gd;
Vp   = sqrt((Ksat + 4/3 * Gsat) ./ rho_sat);
Vs   = sqrt(Gsat ./ rho_sat);
Ip   = (Vp * 1e-3) .* (rho_sat * 1e-3);
VpVs = Vp ./ Vs;
end
%% ---------- Cache curves: exactly identical in CN and EN plots ----------
function V = precompute_visual_curves(A)
V = struct();
V.phi_nodes = 0.05:0.05:0.30;
V.Sg_nodes  = 0:0.25:1;
V.phi_dense = A.cfg.contrast_phi;
V.Sg_dense = A.cfg.contrast_sg;
V.Sat_C = cell(size(V.Sg_nodes)); V.Sat_H = cell(size(V.Sg_nodes));
V.Phi_C = cell(size(V.phi_nodes)); V.Phi_H = cell(size(V.phi_nodes));

Kd_dense = interp1(A.phi_dem_axis, A.K_dem, V.phi_dense, 'pchip');
Gd_dense = interp1(A.phi_dem_axis, A.G_dem, V.phi_dense, 'pchip');
for j = 1:numel(V.Sg_nodes)
    sg = V.Sg_nodes(j);
    [ip, r] = calc_rpt_phi_sg(Kd_dense, Gd_dense, V.phi_dense, ...
        sg, A.T_rpt, A.P_rpt, 'Methane', A.salinity, A.K_solid, A.rho_solid);
    V.Sat_C{j} = [ip(:) r(:)];
    if j > 1 % At Sg=0, the curves coincide with the brine-only curve.
        [ip, r] = calc_rpt_phi_sg(Kd_dense, Gd_dense, V.phi_dense, ...
            sg, A.T_rpt, A.P_rpt, 'Hydrogen', A.salinity, A.K_solid, A.rho_solid);
        V.Sat_H{j} = [ip(:) r(:)];
    end
end
for j = 1:numel(V.phi_nodes)
    phi = V.phi_nodes(j);
    kd = interp1(A.phi_dem_axis, A.K_dem, phi, 'pchip');
    gd = interp1(A.phi_dem_axis, A.G_dem, phi, 'pchip');
    [ip, r] = calc_rpt_phi_sg(kd, gd, phi, V.Sg_dense, ...
        A.T_rpt, A.P_rpt, 'Methane', A.salinity, A.K_solid, A.rho_solid);
    V.Phi_C{j} = [ip(:) r(:)];
    [ip, r] = calc_rpt_phi_sg(kd, gd, phi, V.Sg_dense, ...
        A.T_rpt, A.P_rpt, 'Hydrogen', A.salinity, A.K_solid, A.rho_solid);
    V.Phi_H{j} = [ip(:) r(:)];
end
[V.Ip_c_ref, V.R_c_ref] = calc_rpt_phi_sg(A.Kd_10, A.Gd_10, ...
    A.phi_fixed, .50, A.T_rpt, A.P_rpt, ...
    'Methane', A.salinity, A.K_solid, A.rho_solid);
[V.Ip_h_ref, V.R_h_ref] = calc_rpt_phi_sg(A.Kd_10, A.Gd_10, ...
    A.phi_fixed, .50, A.T_rpt, A.P_rpt, ...
    'Hydrogen', A.salinity, A.K_solid, A.rho_solid);

V.Sg_list = [.20 .50 1.00];
V.P_sweep = A.cfg.contrast_P;
V.T_sweep = A.cfg.contrast_T;
V.C_p = cell(size(V.Sg_list)); V.C_t = cell(size(V.Sg_list));
for j = 1:numel(V.Sg_list)
    sg = V.Sg_list(j);
    [ic, rc] = calc_rpt_phi_sg(A.Kd_10, A.Gd_10, A.phi_fixed, ...
        sg, A.T_rpt, V.P_sweep, 'Methane', ...
        A.salinity, A.K_solid, A.rho_solid);
    [ih, rh] = calc_rpt_phi_sg(A.Kd_10, A.Gd_10, A.phi_fixed, ...
        sg, A.T_rpt, V.P_sweep, 'Hydrogen', ...
        A.salinity, A.K_solid, A.rho_solid);
    V.C_p{j} = calc_ideal_contrast(ih, rh, ic, rc);
    [ic, rc] = calc_rpt_phi_sg(A.Kd_10, A.Gd_10, A.phi_fixed, ...
        sg, V.T_sweep, A.P_rpt, 'Methane', ...
        A.salinity, A.K_solid, A.rho_solid);
    [ih, rh] = calc_rpt_phi_sg(A.Kd_10, A.Gd_10, A.phi_fixed, ...
        sg, V.T_sweep, A.P_rpt, 'Hydrogen', ...
        A.salinity, A.K_solid, A.rho_solid);
    V.C_t{j} = calc_ideal_contrast(ih, rh, ic, rc);
end
end

%% -------------------- Rock physics RPT mesh (render only) ---------------
function [h_sw100, h_ch4, h_h2] = draw_rpt_mesh(ax, V, col_ch4, col_h2, col_sw)
hold(ax, 'on');
h_sw100 = gobjects(1); h_ch4 = gobjects(1); h_h2 = gobjects(1);
for j = 1:numel(V.Sg_nodes)
    C = V.Sat_C{j};
    if j == 1
        h_sw100 = plot(ax, C(:,1), C(:,2), '-', ...
            'Color', col_sw, 'LineWidth', 2.7);
    else
        H = V.Sat_H{j};
        h_ch4 = plot(ax, C(:,1), C(:,2), '--', ...
            'Color', col_ch4, 'LineWidth', 1.55);
        h_h2 = plot(ax, H(:,1), H(:,2), '-', ...
            'Color', col_h2, 'LineWidth', 1.7);
    end
end
for j = 1:numel(V.phi_nodes)
    C = V.Phi_C{j}; H = V.Phi_H{j};
    plot(ax, C(:,1), C(:,2), ':', ...
        'Color', .82*col_ch4 + .18*[1 1 1], ...
        'LineWidth', 1.25, 'HandleVisibility', 'off');
    plot(ax, H(:,1), H(:,2), '-', ...
        'Color', .70*col_h2 + .30*[1 1 1], ...
        'LineWidth', 1.25, 'HandleVisibility', 'off');
end
end

function [Ip, VpVs] = calc_rpt_phi_sg(Kd, Gd, phi, Sg, T_val, P_val, gasType, salinity, K0, rho0)
if ~isnumeric(Sg) || ~isreal(Sg) || isempty(Sg) || ...
        any(~isfinite(Sg(:))) || any(Sg(:) < 0 | Sg(:) > 1)
    error('xu_white_model:InvalidSaturation', 'Gas saturation Sg must be finite and in [0, 1].');
end
[rho_g, K_g] = get_gas_props_vec(gasType, T_val, P_val);
[rho_w, K_w] = get_brine_props(T_val, P_val, salinity);
Sw = 1 - Sg;
Kbrine = gassmann_vec(Kd, K0, phi, K_w);
Kgas   = gassmann_vec(Kd, K0, phi, K_g);
c      = (4/3) .* Gd;
Ksat   = 1 ./ (Sw ./ (Kbrine + c) + Sg ./ (Kgas + c)) - c;
rho_fl = Sg .* rho_g + Sw .* rho_w;
rho_sat = (1 - phi) .* rho0 + phi .* rho_fl;
Gsat    = Gd;
Vp   = sqrt((Ksat + 4/3 * Gsat) ./ rho_sat);
Vs   = sqrt(Gsat ./ rho_sat);
Ip   = (Vp * 1e-3) .* (rho_sat * 1e-3);
VpVs = Vp ./ Vs;
end

function Ksat = gassmann_vec(Kdry, K0, phi, Kfl)
denom = phi ./ Kfl + (1 - phi) ./ K0 - Kdry ./ (K0.^2);
Ksat  = Kdry + ((1 - Kdry ./ K0).^2) ./ denom;
end

function [rho_b, K_b] = get_brine_props(T, P_Pa, S)
P = P_Pa / 1e6;
rho_w = 1 + 1e-6 * (-80 * T - 3.3 * T.^2 + 0.00175 * T.^3 + 489 * P ...
    - 2 * T .* P + 0.016 * T.^2 .* P - 1.3e-5 * T.^3 .* P ...
    - 0.333 * P.^2 - 0.002 * T .* P.^2);
rho_b_gcc = rho_w + S .* (0.668 + 0.44 * S + 1e-6 * (300 * P - 2400 * P .* S ...
            + T .* (80 + 3 * T - 3300 * S - 13 * P + 47 * P .* S)));
w = [1402.85,  4.871, -0.04783,  1.487e-4, -2.197e-7;
     1.524,   -0.0111, 2.747e-4, -6.503e-7,  7.987e-10;
     3.437e-3, 1.739e-4, -2.135e-6, -1.455e-8, 5.230e-11;
    -1.197e-5, -1.628e-6, 1.237e-8, 1.327e-10, -4.614e-13];
Vw = zeros(size(T + P));
for i = 1:5
    for j = 1:4
        Vw = Vw + w(j, i) * (T.^(i-1)) .* (P.^(j-1));
    end
end
Vb = Vw + S .* (1170 - 9.6 * T + 0.055 * T.^2 - 8.5e-5 * T.^3 + 2.6 * P ...
    - 0.0029 * T .* P - 0.0476 * P.^2) ...
     + (S.^1.5) .* (780 - 10 * P + 0.16 * P.^2) - 1820 * (S.^2);
rho_b = rho_b_gcc * 1000;
K_b   = rho_b .* (Vb.^2);
end

function [Kvrh, Gvrh, rhov] = vrh(K, G, rho, V)
Kv   = sum(V .* K);
Kr   = 1 / sum(V ./ K);
Kvrh = (Kv + Kr) / 2;
Gv   = sum(V .* G);
Gr   = 1 / sum(V ./ G);
Gvrh = (Gv + Gr) / 2;
rhov = sum(V .* rho);
end

function [P, Q] = berryman_pq(Km, Gm, Ki, Gi, alpha)
if abs(alpha - 1) < 1e-6
    zeta = Gm / 6 * (9 * Km + 8 * Gm) / (Km + 2 * Gm);
    P = (Km + 4/3 * Gm) / (Ki + 4/3 * Gm);
    Q = (Gm + zeta) / (Gi + zeta);
    return;
end
R = 3 * Gm / (3 * Km + 4 * Gm);
if alpha < 1
    theta = alpha / (1 - alpha^2)^1.5 * (acos(alpha) - alpha * sqrt(1 - alpha^2));
    f     = alpha^2 / (1 - alpha^2) * (3 * theta - 2);
else
    theta = alpha / (alpha^2 - 1)^1.5 * (alpha * sqrt(alpha^2 - 1) - acosh(alpha));
    f     = alpha^2 / (alpha^2 - 1) * (2 - 3 * theta);
end
A  = Gi / Gm - 1;
B  = (1 / 3) * (Ki / Km - Gi / Gm);
F1 = 1 + A * (1.5 * (f + theta) - R * (1.5 * f + 2.5 * theta - 4/3));
F2 = 1 + A * (1 + 1.5 * (f + theta) - 0.5 * R * (3 * f + 5 * theta)) + B * (3 - 4 * R) ...
     + 0.5 * A * (A + 3 * B) * (3 - 4 * R) * (f + theta - R * (f - theta + 2 * theta^2));
F3 = 1 + A * (1 - (f + 1.5 * theta) + R * (f + theta));
F4 = 1 + 0.25 * A * (f + 3 * theta - R * (f - theta));
F5 = A * (-f + R * (f + theta - 4/3)) + B * theta * (3 - 4 * R);
F6 = 1 + A * (1 + f - R * (f + theta)) + B * (1 - theta) * (3 - 4 * R);
F7 = 2 + 0.25 * A * (3 * f + 9 * theta - R * (3 * f + 5 * theta)) + B * theta * (3 - 4 * R);
F8 = A * (1 - 2 * R + 0.5 * f * (R - 1) + 0.5 * theta * (5 * R - 3)) + B * (1 - theta) * (3 - 4 * R);
F9 = A * ((R - 1) * f - R * theta) + B * theta * (3 - 4 * R);
P = F1 / F2;
Q = (1 / 5) * (2 / F3 + 1 / F4 + (F4 * F5 + F6 * F7 - F8 * F9) / (F2 * F4));
end

function [rho, K] = get_gas_props_vec(gasType, T, P)
if ~isnumeric(P) || ~isreal(P) || isempty(P) || ...
        any(~isfinite(P(:))) || any(P(:) <= 0)
    error('xu_white_model:InvalidPressure', 'Absolute pressure P must be finite and greater than zero [Pa].');
end
if ~isnumeric(T) || ~isreal(T) || isempty(T) || ...
        any(~isfinite(T(:))) || any(T(:) <= -273.15)
    error('xu_white_model:InvalidTemperature', 'Temperature T must be finite and above -273.15 degC.');
end
if ~isscalar(T) && ~isscalar(P) && ~isequal(size(T), size(P))
    error('xu_white_model:IncompatibleTPShape', ...
        'T and P must have matching sizes, or one must be scalar. Got T=%s, P=%s.', ...
        mat2str(size(T)), mat2str(size(P)));
end
gasType = validatestring(gasType, {'Methane', 'Hydrogen'}, mfilename, 'gasType');
T_grid = T + zeros(size(P));
P_grid = P + zeros(size(T));
Tk = T_grid + 273.15;
persistent cp_mode
if isempty(cp_mode)
    try
        CoolProp.PropsSI('D', 'T', 333.15, 'P', 10e6, 'Methane');
        CoolProp.PropsSI('A', 'T', 333.15, 'P', 10e6, 'Methane');
        cp_mode = 1;
        fprintf('Gas property backend: MATLAB CoolProp (required).\n');
    catch matlabError
        try
            py.CoolProp.CoolProp.PropsSI('D', 'T', 333.15, 'P', 10e6, 'Methane');
            py.CoolProp.CoolProp.PropsSI('A', 'T', 333.15, 'P', 10e6, 'Methane');
            cp_mode = 2;
            fprintf('Gas property backend: Python CoolProp (required).\n');
        catch pythonError
            dependencyError = MException('xu_white_model:CoolPropUnavailable', ...
                ['CoolProp is required for gas properties. Enable the MATLAB CoolProp wrapper ' ...
                 'or configure pyenv to a Python environment with CoolProp installed. ' ...
                 'No approximate gas equation of state is allowed.']);
            dependencyError = addCause(dependencyError, matlabError);
            dependencyError = addCause(dependencyError, pythonError);
            throwAsCaller(dependencyError);
        end
    end
end
rho = zeros(size(P_grid));
K   = zeros(size(P_grid));
for k = 1:numel(P_grid)
    try
        if cp_mode == 1
            d = double(CoolProp.PropsSI('D', 'T', Tk(k), 'P', P_grid(k), gasType));
            c = double(CoolProp.PropsSI('A', 'T', Tk(k), 'P', P_grid(k), gasType));
        else
            d = double(py.CoolProp.CoolProp.PropsSI('D', 'T', Tk(k), 'P', P_grid(k), gasType));
            c = double(py.CoolProp.CoolProp.PropsSI('A', 'T', Tk(k), 'P', P_grid(k), gasType));
        end
    catch propertyError
        stateError = MException('xu_white_model:CoolPropStateFailure', ...
            'CoolProp failed for %s at T=%.8g degC, P=%.8g MPa (absolute).', ...
            gasType, T_grid(k), P_grid(k) / 1e6);
        throwAsCaller(addCause(stateError, propertyError));
    end
    bulk = d * c^2;
    if ~isreal(d) || ~isreal(c) || ~isfinite(d) || ~isfinite(c) || ...
            d <= 0 || c <= 0 || ~isfinite(bulk) || bulk <= 0
        error('xu_white_model:CoolPropInvalidProperties', ...
            'Invalid CoolProp result for %s at T=%.8g degC, P=%.8g MPa: rho=%g, c=%g, K=%g.', ...
            gasType, T_grid(k), P_grid(k) / 1e6, d, c, bulk);
    end
    rho(k) = d;
    K(k)   = bulk;
end
end


function validate_test3_settings(cfg)
% Reject silently clipped T/P windows; both A and B have fixed true support.
assert(any(strcmp(cfg.normal_sampling, {'clipped','truncated'})), ...
    'RPT:SamplingModel', 'normal_sampling must be clipped or truncated.');
assert(cfg.n_map >= 16 && cfg.n_map == fix(cfg.n_map) && ...
    cfg.n_reference >= cfg.n_map && cfg.n_reference == fix(cfg.n_reference), ...
    'RPT:SampleCount', 'Use integer sample counts: 16 <= n_map <= n_reference.');
assert(cfg.batch_pixels >= 1 && cfg.batch_pixels == fix(cfg.batch_pixels), ...
    'RPT:BatchSize', 'batch_pixels must be a positive integer.');
assert(cfg.phi_fixed > 0 && cfg.phi_fixed <= .35 && ...
    cfg.sg_fixed >= 0 && cfg.sg_fixed <= 1 && ...
    cfg.sigma_phi >= 0 && cfg.sigma_sg >= 0 && ...
    cfg.delta_T >= 0 && cfg.delta_P >= 0, ...
    'RPT:PhysicalParameters', 'Invalid fixed state or nuisance width.');
axes_fields = {'A_phi','A_sg','A_T','A_P','B_sigma_phi','B_sigma_sg', ...
    'B_delta_T','B_delta_P','contrast_phi','contrast_sg','contrast_T','contrast_P', ...
    'lut_T','lut_P'};
for k = 1:numel(axes_fields)
    v = cfg.(axes_fields{k});
    assert(isnumeric(v) && isreal(v) && isvector(v) && numel(v) >= 2 && ...
        all(isfinite(v)) && all(diff(v) > 0), 'RPT:ParameterAxis', ...
        '%s must contain at least two increasing finite coordinates.', axes_fields{k});
end
assert(all(cfg.A_phi >= cfg.phi_bounds(1) & cfg.A_phi <= cfg.phi_bounds(2)) && ...
    all(cfg.A_sg >= cfg.sg_bounds(1) & cfg.A_sg <= cfg.sg_bounds(2)) && ...
    all(cfg.B_sigma_phi >= 0) && all(cfg.B_sigma_sg >= 0) && ...
    all(cfg.B_delta_T >= 0) && all(cfg.B_delta_P >= 0), ...
    'RPT:EnsembleDomain', 'Centers/widths lie outside the configured domain.');
tp_extremes = [min(cfg.A_T)-cfg.delta_T max(cfg.A_T)+cfg.delta_T; ...
    cfg.T_center-max(cfg.B_delta_T) cfg.T_center+max(cfg.B_delta_T)];
pp_extremes = [min(cfg.A_P)-cfg.delta_P max(cfg.A_P)+cfg.delta_P; ...
    cfg.P_center-max(cfg.B_delta_P) cfg.P_center+max(cfg.B_delta_P)];
tol_T = 1e-10; tol_P = 1e-4;
assert(all(tp_extremes(:,1) >= cfg.T_bounds(1)-tol_T) && ...
    all(tp_extremes(:,2) <= cfg.T_bounds(2)+tol_T) && ...
    all(pp_extremes(:,1) >= cfg.P_bounds(1)-tol_P) && ...
    all(pp_extremes(:,2) <= cfg.P_bounds(2)+tol_P), ...
    'RPT:WindowOutOfRange', ...
    'A/B windows exceed T/P support. Adjust centers or widths explicitly; no window shrinking.');
assert(min(cfg.lut_T) <= min([cfg.T_bounds cfg.contrast_T cfg.T_ref]) && ...
    max(cfg.lut_T) >= max([cfg.T_bounds cfg.contrast_T cfg.T_ref]) && ...
    min(cfg.lut_P) <= min([cfg.P_bounds cfg.contrast_P cfg.P_ref]) && ...
    max(cfg.lut_P) >= max([cfg.P_bounds cfg.contrast_P cfg.P_ref]), ...
    'RPT:LUTCoverage', 'CoolProp LUT must cover all fixed states and windows.');
assert(iscell(cfg.languages) && ~isempty(cfg.languages) && ...
    all(ismember(cfg.languages, {'CN','EN'})), ...
    'RPT:Language', 'languages must contain CN and/or EN.');
end

function U = shared_sample_design(n, seed)
% Common shifted Halton points; nested prefixes permit N->2N comparisons.
% No Statistics Toolbox; preserve caller RNG state. This is not added noise.
assert(n >= 1 && n == fix(n), 'RPT:SampleCount', 'n must be a positive integer.');
old_rng = rng;
rng(seed, 'twister'); shift = rand(1,2);
rng(old_rng);
U = zeros(n,2);
for d = 1:2
    base = d + 1; indices = (1:n)'; factor = 1/base;
    while any(indices)
        U(:,d) = U(:,d) + factor*mod(indices,base);
        indices = floor(indices/base); factor = factor/base;
    end
end
U = mod(bsxfun(@plus, U, shift), 1);
U = min(max(U, eps), 1-eps);
end

function [phi, Sg, T, P] = nuisance_samples(kind, points, cfg, U)
% Columns are pixels, rows are shared samples. P and deltaP are in Pa.
assert(size(points,2) == 2 && all(isfinite(points(:))) && ...
    size(U,2) == 2 && all(U(:) > 0 & U(:) < 1), ...
    'RPT:EnsembleInput', 'Expected finite M-by-2 points and open-unit-square samples.');
n = size(U,1); m = size(points,1); one = ones(n,m);
switch kind
    case 'A_phi'
        phi = bounded_normal_samples(U(:,1), points(:,1)', ...
            cfg.sigma_phi, cfg.phi_bounds, cfg.normal_sampling);
        Sg = bounded_normal_samples(U(:,2), points(:,2)', ...
            cfg.sigma_sg, cfg.sg_bounds, cfg.normal_sampling);
        T = cfg.T_ref*one; P = cfg.P_ref*one;
    case 'B_phi'
        phi = bounded_normal_samples(U(:,1), cfg.phi_ref, ...
            points(:,1)', cfg.phi_bounds, cfg.normal_sampling);
        Sg = bounded_normal_samples(U(:,2), cfg.sg_ref, ...
            points(:,2)', cfg.sg_bounds, cfg.normal_sampling);
        T = cfg.T_ref*one; P = cfg.P_ref*one;
    case 'A_tp'
        T = bsxfun(@plus, cfg.delta_T*(2*U(:,1)-1), points(:,2)');
        P = bsxfun(@plus, cfg.delta_P*(2*U(:,2)-1), points(:,1)');
        phi = cfg.phi_fixed*one; Sg = cfg.sg_fixed*one;
    case 'B_tp'
        T = cfg.T_center + bsxfun(@times, 2*U(:,1)-1, points(:,2)');
        P = cfg.P_center + bsxfun(@times, 2*U(:,2)-1, points(:,1)');
        phi = cfg.phi_fixed*one; Sg = cfg.sg_fixed*one;
    otherwise
        error('RPT:UnknownMap', 'Unknown ensemble kind %s.', kind);
end
assert(all(T(:) >= cfg.T_bounds(1)-1e-10 & T(:) <= cfg.T_bounds(2)+1e-10) && ...
    all(P(:) >= cfg.P_bounds(1)-1e-4 & P(:) <= cfg.P_bounds(2)+1e-4), ...
    'RPT:SampleOutOfRange', 'T/P samples exceed support; never silently truncate uniform windows.');
end

function X = bounded_normal_samples(u, mu, sigma, bounds, sampling)
% mu/sigma are nominal, BEFORE clipping/truncation; output moments differ.
mu = mu(:)'; sigma = sigma(:)';
m = max(numel(mu),numel(sigma));
assert((isscalar(mu) || numel(mu)==m) && (isscalar(sigma) || numel(sigma)==m) && ...
    all(isfinite(mu)) && all(isfinite(sigma)) && all(sigma >= 0) && ...
    all(mu >= bounds(1) & mu <= bounds(2)), ...
    'RPT:NormalParameters', 'Invalid normal center, width, or support.');
mu = mu + zeros(1,m); sigma = sigma + zeros(1,m);
if strcmp(sampling, 'clipped')
    z = -sqrt(2)*erfcinv(2*u(:));
    X = bsxfun(@plus, bsxfun(@times, z, sigma), mu);
    X = min(max(X, bounds(1)), bounds(2));
elseif strcmp(sampling, 'truncated')
    X = repmat(mu,numel(u),1);
    active = sigma > 0;
    a = .5*erfc(-(bounds(1)-mu(active))./sigma(active)/sqrt(2));
    b = .5*erfc(-(bounds(2)-mu(active))./sigma(active)/sqrt(2));
    q = bsxfun(@plus, bsxfun(@times, u(:), b-a), a);
    z = -sqrt(2)*erfcinv(2*min(max(q,eps),1-eps));
    X(:,active) = bsxfun(@plus, bsxfun(@times,z,sigma(active)),mu(active));
else
    error('RPT:SamplingModel', 'Unknown normal sampling model: %s.', sampling);
end
end

function [values, clouds] = eval_separability_at(kind, points, A, lut, cfg, U, opt)
% Vectorized paired forward calculations, bounded memory, no parfor needed.
if nargin < 7 || isempty(opt), opt = cfg.metric; end
count = size(points,1);
names = {'S_ellipse','S_lda','S_gaussian2d','S_kde','Pe_gaussian2d','delta', ...
    'actual_sigma_phi','actual_sigma_sg','integration_error','rank_c','rank_h'};
values = struct();
for k = 1:numel(names), values.(names{k}) = nan(count,1); end
if nargout > 1, clouds = cell(count,2); end
for first = 1:cfg.batch_pixels:count
    rows = first:min(first+cfg.batch_pixels-1,count);
    [phi,Sg,T,P] = nuisance_samples(kind,points(rows,:),cfg,U);
    Kd = interp1(A.phi_dem_axis,A.K_dem,phi,'pchip');
    Gd = interp1(A.phi_dem_axis,A.G_dem,phi,'pchip');
    [Ic,Rc] = calc_rpt_lut(Kd,Gd,phi,Sg,T,P,'Methane', ...
        A.salinity,A.K_solid,A.rho_solid,lut);
    [Ih,Rh] = calc_rpt_lut(Kd,Gd,phi,Sg,T,P,'Hydrogen', ...
        A.salinity,A.K_solid,A.rho_solid,lut);
    for j = 1:numel(rows)
        index = rows(j); Xc = [Ic(:,j) Rc(:,j)]; Xh = [Ih(:,j) Rh(:,j)];
        st = ensemble_separability(Xc,Xh,opt);
        for f = 1:6, values.(names{f})(index) = st.(names{f}); end
        values.rank_c(index) = st.rank_c; values.rank_h(index) = st.rank_h;
        values.actual_sigma_phi(index) = std(phi(:,j));
        values.actual_sigma_sg(index) = std(Sg(:,j));
        values.integration_error(index) = st.diagnostics.gaussian2d_error;
        if nargout > 1, clouds(index,:) = {Xc,Xh}; end
    end
end
end

function map = build_separability_map(kind,A,lut,cfg,U)
map = separability_map_metadata(kind,cfg);
[X,Y] = meshgrid(map.x_values,map.y_values);
points = [X(:) Y(:)]; count = size(points,1);
fprintf('\n%s: %d x %d pixels, %d paired samples/pixel.\n', ...
    map.name,numel(map.x),numel(map.y),size(U,1));
started = tic;
names = {'S_ellipse','S_lda','S_gaussian2d','S_kde','Pe_gaussian2d','delta', ...
    'actual_sigma_phi','actual_sigma_sg','integration_error','rank_c','rank_h'};
for k = 1:numel(names), map.(names{k}) = nan(size(X)); end
% Progress is printed about20 times per map; intermediate original figures
% remain responsive. Block size changes neither samples nor numerical scores.
chunk = max(cfg.batch_pixels,ceil(count/20));
for first = 1:chunk:count
    rows = first:min(first+chunk-1,count);
    values = eval_separability_at(kind,points(rows,:),A,lut,cfg,U,cfg.metric);
    for k = 1:numel(names), map.(names{k})(rows) = values.(names{k}); end
    fprintf('  %s: %5.1f%%, %.1f s elapsed\n',kind,100*rows(end)/count,toc(started));
    drawnow;
end
map.samples = size(U,1);
map.elapsed_seconds = toc(started);
map.normal_sampling = cfg.normal_sampling;
map.near_complete_threshold = .99;
map.Pe_lda = (1-map.S_lda)/2;
map.Pe_kde = (1-map.S_kde)/2;
% Retain these diagnostic arrays for manual alternative color scales. A
% value of Inf is an exact/underflow-zero error, never an invented finite log.
map.minus_log10_Pe = -log10(map.Pe_gaussian2d);
end

function map = separability_map_metadata(kind,cfg)
map = struct('kind',kind,'reference_xy',[]);
% Explicit labels prevent nominal sigma from being read as realized sigma.
if strcmp(cfg.normal_sampling,'clipped')
    sampling_cn = '截边正态；σ 为截边前宽度';
    sampling_en = 'Clipped normal; sigma is the pre-clipping width';
else
    sampling_cn = '截断正态；σ 为截断前宽度';
    sampling_en = 'Truncated normal; sigma is the pre-truncation width';
end
switch kind
    case 'A_phi'
        map.name = 'Figure 3 - A: phi/Sg centers';
        map.x_values = cfg.A_phi; map.y_values = cfg.A_sg;
        map.x = 100*map.x_values; map.y = 100*map.y_values;
        map.xlabel_cn = '孔隙度中心 φ (%)'; map.xlabel_en = 'Porosity center φ (%)';
        map.ylabel_cn = '气体饱和度中心 S_g (%)'; map.ylabel_en = 'Gas saturation center S_g (%)';
        common = sprintf('T = %g °C, P = %g MPa; σφ = %g pp, σSg = %g pp', ...
            cfg.T_ref,cfg.P_ref/1e6,100*cfg.sigma_phi,100*cfg.sigma_sg);
        map.note_cn = {common; [sampling_cn '；不添加观测噪声']};
        map.note_en = {common; [sampling_en '; no observation noise']};
    case 'A_tp'
        map.name = 'Figure 4 - A: T/P centers';
        map.x_values = cfg.A_P; map.y_values = cfg.A_T;
        map.x = map.x_values/1e6; map.y = map.y_values;
        map.xlabel_cn = '压力中心 P (MPa)'; map.xlabel_en = 'Pressure center P (MPa)';
        map.ylabel_cn = '温度中心 T (°C)'; map.ylabel_en = 'Temperature center T (°C)';
        common = sprintf('φ = %g%%, S_g = %g%%; ΔT = %g °C, ΔP = %g MPa', ...
            100*cfg.phi_fixed,100*cfg.sg_fixed,cfg.delta_T,cfg.delta_P/1e6);
        map.note_cn = {common; '中心 ± 固定半宽的均匀分布；不添加观测噪声'};
        map.note_en = {common; 'Uniform center +/- fixed half widths; no observation noise'};
    case 'B_phi'
        map.name = 'Figure 5 - B: phi/Sg widths';
        map.x_values = cfg.B_sigma_phi; map.y_values = cfg.B_sigma_sg;
        map.x = 100*map.x_values; map.y = 100*map.y_values;
        map.xlabel_cn = '孔隙度标准差 σφ (百分点)'; map.xlabel_en = 'Porosity sigma (percentage points)';
        map.ylabel_cn = '饱和度标准差 σSg (百分点)'; map.ylabel_en = 'Saturation sigma (percentage points)';
        common = sprintf('φ = %g%%, S_g = %g%%; T = %g °C, P = %g MPa', ...
            100*cfg.phi_ref,100*cfg.sg_ref,cfg.T_ref,cfg.P_ref/1e6);
        map.note_cn = {common; [sampling_cn '；参考参数来自 Fig2 上排']};
        map.note_en = {common; [sampling_en '; reference: Fig2 top-row widths']};
        map.reference_xy = 100*[cfg.sigma_phi cfg.sigma_sg];
    case 'B_tp'
        map.name = 'Figure 6 - B: T/P widths';
        map.x_values = cfg.B_delta_P; map.y_values = cfg.B_delta_T;
        map.x = map.x_values/1e6; map.y = map.y_values;
        map.xlabel_cn = '压力半宽 ΔP (MPa)'; map.xlabel_en = 'Pressure half width ΔP (MPa)';
        map.ylabel_cn = '温度半宽 ΔT (°C)'; map.ylabel_en = 'Temperature half width ΔT (°C)';
        common = sprintf('φ = %g%%, S_g = %g%%; T_c = %g °C, P_c = %g MPa', ...
            100*cfg.phi_fixed,100*cfg.sg_fixed,cfg.T_center,cfg.P_center/1e6);
        map.note_cn = {common; '均匀分布半宽；参考参数来自 Fig2 下排的全域窗口'};
        map.note_en = {common; 'Uniform half widths; reference: Fig2 bottom-row full window'};
        if abs(cfg.T_center-mean(cfg.T_bounds))<1e-10 && ...
                abs(cfg.P_center-mean(cfg.P_bounds))<1e-4
            map.reference_xy = [diff(cfg.P_bounds)/2e6 diff(cfg.T_bounds)/2];
        else
            map.note_cn{2} = '均匀分布半宽；中心与 Fig2 不同，不显示参考星号';
            map.note_en{2} = 'Uniform half widths; center differs from Fig2, so no reference star';
        end
    otherwise
        error('RPT:UnknownMap', 'Unknown map kind %s.',kind);
end
map.x_values = map.x_values(:)'; map.y_values = map.y_values(:)';
if ~isempty(map.reference_xy) && ...
        (map.reference_xy(1)<min(map.x) || map.reference_xy(1)>max(map.x) || ...
         map.reference_xy(2)<min(map.y) || map.reference_xy(2)>max(map.y))
    map.reference_xy = [];
    map.note_cn{end} = [map.note_cn{end} '；参考点不在当前坐标范围内'];
    map.note_en{end} = [map.note_en{end} '; reference is outside the displayed axes'];
end
map.reference_label_cn = 'Fig2 参考'; map.reference_label_en = 'Fig2 reference';
end

function report = check_separability_resolution(map,A,lut,cfg)
% Probe actual NEW midpoints, not unchanged shared nodes. This is a local
% convergence diagnostic, not a proof over every pixel or a smoothing step.
fields = {'S_ellipse','S_lda','S_gaussian2d','S_kde'};
nx = numel(map.x_values); ny = numel(map.y_values);
cells = [max(1,floor(ny/2)) max(1,floor(nx/2)); ...
    max(1,floor(ny/4)) max(1,floor(nx/4)); ...
    min(ny-1,max(1,floor(3*ny/4))) min(nx-1,max(1,floor(3*nx/4)))];
% Add cells closest to important probability contours; ignore invalid cells.
for f = [2 3]
    V = map.(fields{f});
    mid = (V(1:end-1,1:end-1)+V(2:end,1:end-1)+ ...
        V(1:end-1,2:end)+V(2:end,2:end))/4;
    for level = [.90 .99]
        distance = abs(mid-level); distance(~isfinite(distance)) = Inf;
        [nearest,index] = min(distance(:));
        if isfinite(nearest)
            [i,j] = ind2sub(size(mid),index); cells(end+1,:) = [i j]; %#ok<AGROW>
        end
    end
end
cells = unique(cells,'rows');
points = [(map.x_values(cells(:,2))+map.x_values(cells(:,2)+1))'/2, ...
    (map.y_values(cells(:,1))+map.y_values(cells(:,1)+1))'/2];
% Axis row/column shape must not change the number of probe points.
points = reshape(points,[],2);
U = shared_sample_design(cfg.n_map,cfg.seed);
base = eval_separability_at(map.kind,points,A,lut,cfg,U,cfg.metric);
fine_samples = eval_separability_at(map.kind,points,A,lut,cfg, ...
    shared_sample_design(2*cfg.n_map,cfg.seed),cfg.metric);
fine_opt = cfg.metric;
fine_opt.ellipse_points = 2*cfg.metric.ellipse_points;
fine_opt.kde_points = 2*cfg.metric.kde_points;
fine_opt.quadrature_order = 2*cfg.metric.quadrature_order;
fine_method = eval_separability_at(map.kind,points,A,lut,cfg,U,fine_opt);
report = struct('points',points,'methods',{fields}, ...
    'sample_error',nan(1,4),'method_error',nan(1,4),'midpoint_error',nan(1,4), ...
    'tolerance',cfg.convergence_tolerance,'scope','Selected midpoint/contour probes only');
for f = 1:numel(fields)
    name = fields{f};
    interpolated = interp2(map.x_values,map.y_values,map.(name), ...
        points(:,1),points(:,2),'linear');
    report.sample_error(f) = finite_max(abs(base.(name)-fine_samples.(name)));
    report.method_error(f) = finite_max(abs(base.(name)-fine_method.(name)));
    report.midpoint_error(f) = finite_max(abs(base.(name)-interpolated));
end
fprintf('  %s probe errors [ellipse LDA-G 2D-G KDE], percentage points:\n',map.kind);
fprintf('    sample doubling: %s\n',mat2str(100*report.sample_error,3));
fprintf('    method refinement: %s\n',mat2str(100*report.method_error,3));
fprintf('    true midpoint vs grid interpolation: %s\n',mat2str(100*report.midpoint_error,3));
worst = finite_max([report.sample_error report.method_error report.midpoint_error]);
report.exceeds_tolerance = isfinite(worst) && worst > cfg.convergence_tolerance;
if report.exceeds_tolerance
    warning('RPT:ResolutionCheck', ...
        '%s has a probe change > %.2f percentage points. Increase n_map/method resolution or refine parameter steps.', ...
        map.kind,100*cfg.convergence_tolerance);
end
end

function report = check_lut_references(A,lut,cfg,U)
% Compare identical physical samples via LUT and DIRECT CoolProp. This
% checks numerical interpolation, without fitting any observation noise.
U = U(1:min(256,size(U,1)),:);
kinds = {'B_phi','B_tp'};
points = {[cfg.sigma_phi cfg.sigma_sg], ...
    [diff(cfg.P_bounds)/2 diff(cfg.T_bounds)/2]};
report = struct('kind',kinds,'max_relative_Ip',[],'max_relative_ratio',[], ...
    'score_difference',[]);
for k = 1:2
    local_cfg = cfg;
    if k==2, local_cfg = reference_tp_settings(cfg); end
    [phi,Sg,T,P] = nuisance_samples(kinds{k},points{k},local_cfg,U);
    Kd = interp1(A.phi_dem_axis,A.K_dem,phi,'pchip');
    Gd = interp1(A.phi_dem_axis,A.G_dem,phi,'pchip');
    [ic,rc] = calc_rpt_lut(Kd,Gd,phi,Sg,T,P,'Methane',A.salinity,A.K_solid,A.rho_solid,lut);
    [ih,rh] = calc_rpt_lut(Kd,Gd,phi,Sg,T,P,'Hydrogen',A.salinity,A.K_solid,A.rho_solid,lut);
    [dc,ec] = calc_rpt_phi_sg(Kd,Gd,phi,Sg,T,P,'Methane',A.salinity,A.K_solid,A.rho_solid);
    [dh,eh] = calc_rpt_phi_sg(Kd,Gd,phi,Sg,T,P,'Hydrogen',A.salinity,A.K_solid,A.rho_solid);
    sl = ensemble_separability([ic rc],[ih rh],cfg.metric);
    sd = ensemble_separability([dc ec],[dh eh],cfg.metric);
    report(k).max_relative_Ip = max(abs([ic-dc;ih-dh])./max(abs([dc;dh]),realmin));
    report(k).max_relative_ratio = max(abs([rc-ec;rh-eh])./max(abs([ec;eh]),realmin));
    report(k).score_difference = abs([sl.S_ellipse sl.S_lda sl.S_gaussian2d sl.S_kde] - ...
        [sd.S_ellipse sd.S_lda sd.S_gaussian2d sd.S_kde]);
    fprintf('LUT check %s: max relative Ip %.3g, ratio %.3g; score changes(pp) %s\n', ...
        kinds{k},report(k).max_relative_Ip,report(k).max_relative_ratio, ...
        mat2str(100*report(k).score_difference,3));
    if finite_max(report(k).score_difference) > cfg.convergence_tolerance
        warning('RPT:LUTResolution', ...
            '%s LUT vs direct CoolProp exceeds score tolerance; refine cfg.lut_T/lut_P.',kinds{k});
    end
end
end

function print_reference_scores(name,st)
fprintf(['%s: ellipse-area %.3f%%; LDA Gaussian %.3f%%; 2D Gaussian %.3f%%; LDA KDE %.3f%%.\n' ...
    '  Equal-prior Pe: LDA-G %.6g%%; 2D-G %.6g%%; KDE %.6g%% (ellipse has no Pe).\n'], ...
    name,100*st.S_ellipse,100*st.S_lda,100*st.S_gaussian2d,100*st.S_kde, ...
    100*st.Pe_lda,100*st.Pe_gaussian2d,100*st.Pe_kde);
end

function value = finite_max(values)
values = values(isfinite(values));
if isempty(values), value = NaN; else, value = max(values(:)); end
end


function [cn,en] = fig2_condition_text(cfg,row)
if row == 1
    cn = sprintf('固定 T = %g °C、P = %g MPa；仅 φ、S_g 变化',cfg.T_ref,cfg.P_ref/1e6);
    en = sprintf('Fixed T = %g °C, P = %g MPa; variable φ and S_g',cfg.T_ref,cfg.P_ref/1e6);
else
    cn = sprintf('固定 φ = %g%%、S_g = %g%%；仅 T、P 变化',100*cfg.phi_fixed,100*cfg.sg_fixed);
    en = sprintf('Fixed φ = %g%%, S_g = %g%%; variable T and P',100*cfg.phi_fixed,100*cfg.sg_fixed);
end
end

function ref = reference_tp_settings(cfg)
% Fig2 always spans the specified physical support, independently of the
% freely adjustable center of the B tolerance map.
ref = cfg;
ref.T_center = mean(cfg.T_bounds); ref.P_center = mean(cfg.P_bounds);
end

function [grid,fc,fh,atoms] = reference_projection_densities(st)
% Reuse the SAME KDE curves as the reported Pe. No second coarse integration.
% Empty cached curves arise for atoms/identical/disjoint KDEs; then create
% display-only curves on class-resolved grids, keeping st.Pe_kde unchanged.
c = st.z_c(:); h = st.z_h(:);
if isempty(c) || isempty(h)
    grid=[]; fc=[]; fh=[]; atoms=[false false]; return;
end
atoms = [max(c)==min(c), max(h)==min(h)];
if ~isempty(st.kde_grid)
    grid=st.kde_grid(:); fc=st.kde_pdf_c(:); fh=st.kde_pdf_h(:); return;
end
samples = {c,h}; parts=cell(1,2);
for k=1:2
    z=samples{k};
    if atoms(k)
        parts{k}=z(1);
    else
        bw=1.06*std(z)*numel(z)^(-.2);
        parts{k}=linspace(min(z)-8*bw,max(z)+8*bw,512)';
    end
end
grid=unique([parts{1}(:);parts{2}(:)]);
fc=zeros(size(grid)); fh=zeros(size(grid));
if ~atoms(1), fc=calc_1d_kde(c,grid)'; end
if ~atoms(2), fh=calc_1d_kde(h,grid)'; end
end

function draw_reference_marginal(ax,c,h,grid,col_c,col_h,orientation)
% Continuous densities retain the original KDE; an exact point mass is a
% line marker, never a tiny artificial Gaussian bandwidth or a zero ylim.
samples={c(:),h(:)}; colors={col_c,col_h}; peak=0;
hold(ax,'on');
for k=1:2
    z=samples{k}; color=colors{k};
    if max(z)==min(z)
        if strcmp(orientation,'top')
            xline(ax,z(1),'-','Color',color,'LineWidth',2,'HandleVisibility','off');
        else
            yline(ax,z(1),'-','Color',color,'LineWidth',2,'HandleVisibility','off');
        end
    else
        density=calc_1d_kde(z,grid); peak=max(peak,max(density));
        if strcmp(orientation,'top')
            fill(ax,grid,density,.24*color+.76*[1 1 1], ...
                'EdgeColor',color,'LineWidth',1.9,'HandleVisibility','off');
        else
            fill(ax,density,grid,.24*color+.76*[1 1 1], ...
                'EdgeColor',color,'LineWidth',1.9,'HandleVisibility','off');
        end
    end
end
if peak<=0, peak=1; end
if strcmp(orientation,'top')
    xlim(ax,[grid(1) grid(end)]); ylim(ax,[0 1.06*peak]);
else
    xlim(ax,[0 1.06*peak]); ylim(ax,[grid(1) grid(end)]);
end
axis(ax,'off');
end


function stats = ensemble_separability(Xc, Xh, opt)
% Four deliberately different separation measures for two N-by-2 ensembles.
% No observational noise or dimensional covariance ridge is introduced.
% Covariances use N-1, matching MATLAB cov(X). Both classes share feature units.
% S_ellipse is 1 - intersection/union of 95% covariance ellipses (a heuristic).
% S_lda and S_gaussian2d are total-variation distances of the stated Gaussian
% models. Pe=(1-S)/2 assumes equal class priors. S_kde uses the SAME projection.
% Distinct singular supports are mutually singular under exact observation.
if nargin<3, opt=struct(); end
opt=separability_options(opt);
assert(isnumeric(Xc) && isnumeric(Xh) && isreal(Xc) && isreal(Xh) && ...
    size(Xc,2)==2 && size(Xh,2)==2 && ~isempty(Xc) && ~isempty(Xh) && ...
    all(isfinite(Xc(:))) && all(isfinite(Xh(:))), ...
    'RPT:InvalidEnsemble','Each class must be a finite real N-by-2 array.');

% Center before averaging: constant clouds must have exactly zero covariance.
Dc0=bsxfun(@minus,Xc,Xc(1,:)); Dh0=bsxfun(@minus,Xh,Xh(1,:));
ac=mean(Dc0,1); ah=mean(Dh0,1);
Dc=bsxfun(@minus,Dc0,ac); Dh=bsxfun(@minus,Dh0,ah);
mc=Xc(1,:)+ac; mh=Xh(1,:)+ah;
d=(Xh(1,:)-Xc(1,:))+(ah-ac);
Cc=(Dc'*Dc)/max(size(Xc,1)-1,1);
Ch=(Dh'*Dh)/max(size(Xh,1)-1,1);
units=max([sqrt(max(diag(Cc),0))';sqrt(max(diag(Ch),0))';abs(d)],[],1);
units(units==0)=1;
Dc=bsxfun(@rdivide,Dc,units); Dh=bsxfun(@rdivide,Dh,units);
Ccs=(Dc'*Dc)/max(size(Xc,1)-1,1);
Chs=(Dh'*Dh)/max(size(Xh,1)-1,1);
ds=d./units;
[Vc,lc,rc]=covariance_axes(Ccs,opt.rank_tol);
[Vh,lh,rh]=covariance_axes(Chs,opt.rank_tol);
Cp=.5*(Ccs+Chs);
[Vp,lp,rp]=covariance_axes(Cp,opt.rank_tol);
stats=struct('S_ellipse',NaN,'S_lda',NaN,'Pe_lda',NaN, ...
    'S_gaussian2d',NaN,'Pe_gaussian2d',NaN,'S_kde',NaN,'Pe_kde',NaN, ...
    'delta',NaN,'z_c',[],'z_h',[],'direction',[NaN;NaN], ...
    'mu_c',mc,'mu_h',mh,'cov_c',Cc,'cov_h',Ch, ...
    'rank_c',rc,'rank_h',rh,'rank_pooled',rp, ...
    'm_c',NaN,'m_h',NaN,'s_c',NaN,'s_h',NaN, ...
    'kde_grid',[],'kde_pdf_c',[],'kde_pdf_h',[]);
diagout=struct('lda_defined',true,'projection_mode','Fisher', ...
    'feature_scale',units,'pooled_eigenvalues',lp', ...
    'gaussian2d_error',0,'gaussian2d_tail_bound',0, ...
    'gaussian2d_method','','kde_points',0,'kde_tail_bound',NaN, ...
    'kde_grid_step',NaN,'kde_bandwidths',[NaN NaN], ...
    'kde_resolved',false,'ellipse_defined',rc==2 && rh==2, ...
    'ellipse_points',opt.ellipse_points);

same_mean=norm(ds)<=opt.support_tol;
same_cov=norm(Ccs-Chs,'fro')<=opt.rank_tol*max(1,norm(Cp,'fro'));
null_shift=Vp(:,lp==0)'*ds';
if same_mean && ~same_cov
    % A variance-only discrimination problem has no Fisher mean direction.
    % Leave the projected methods undefined; 2D Gaussian remains meaningful.
    diagout.lda_defined=false;
    diagout.projection_mode='undefined: equal means, different covariances';
    stats.delta=0;
else
    null_separation=~same_mean && norm(null_shift)>opt.support_tol*(1+norm(ds));
    if null_separation
        w=Vp(:,lp==0)*null_shift;
        diagout.projection_mode='Fisher null-space separation';
    elseif same_mean
        w=Vp(:,end);
        diagout.lda_defined=false;
        diagout.projection_mode='identical Gaussian moments: fixed principal axis';
    else
        w=Vp(:,lp>0)*((Vp(:,lp>0)'*ds')./lp(lp>0));
    end
    if norm(w)==0, w=[1;0]; end
    w=w/norm(w);
    if ds*w<0, w=-w; end
    sep=ds*w;
    pc=Dc*w; ph=Dh*w;
    sc=sqrt(sum(pc.^2)/max(size(Xc,1)-1,1));
    sh=sqrt(sum(ph.^2)/max(size(Xh,1)-1,1));
    if null_separation, sc=0; sh=0; pc(:)=0; ph(:)=0; end
    ps=sqrt(.5*(sc^2+sh^2));
    if ps==0
        pscale=max(abs(sep),1);
        if same_mean, stats.delta=0; else, stats.delta=Inf; end
    else
        pscale=ps; stats.delta=abs(sep)/ps;
    end
    stats.direction=bsxfun(@rdivide,w,units')/pscale;
    stats.z_c=(pc-.5*sep)/pscale; stats.z_h=(ph+.5*sep)/pscale;
    stats.m_c=-.5*sep/pscale; stats.m_h=.5*sep/pscale;
    stats.s_c=sc/pscale; stats.s_h=sh/pscale;
    ov=gaussian_overlap_1d(stats.m_c,stats.s_c,stats.m_h,stats.s_h);
    stats.S_lda=1-ov; stats.Pe_lda=.5*ov;
    % Identical Gaussian moments leave no distinguished projection for different
    % empirical shapes. Exact identical clouds/atoms still have overlap one.
    if opt.compute_kde && (~same_mean || isequal(Xc,Xh) || (sc==0 && sh==0))
        [ov,grid,fc,fh,kdiag]=projected_kde_overlap( ...
            stats.z_c,stats.z_h,stats.s_c,stats.s_h,opt);
        stats.S_kde=1-ov; stats.Pe_kde=.5*ov;
        stats.kde_grid=grid; stats.kde_pdf_c=fc; stats.kde_pdf_h=fh;
        diagout.kde_points=kdiag.points; diagout.kde_grid_step=kdiag.step;
        diagout.kde_bandwidths=kdiag.bandwidths;
        diagout.kde_tail_bound=kdiag.tail_bound;
        diagout.kde_resolved=kdiag.resolved;
    end
end

[ov,err,tail,method]=gaussian_overlap_2d(ds,Ccs,Chs,Vc,lc,rc,Vh,lh,rh,Vp,lp,opt);
stats.S_gaussian2d=1-ov; stats.Pe_gaussian2d=.5*ov;
diagout.gaussian2d_error=err; diagout.gaussian2d_tail_bound=tail;
diagout.gaussian2d_method=method;
if rc==2 && rh==2
    stats.S_ellipse=ellipse_jaccard_separation(-ds/2,Ccs,ds/2,Chs,opt.ellipse_points);
end
stats.diagnostics=diagout;
end

function opt=separability_options(opt)
defaults=struct('ellipse_points',128,'quadrature_order',24,'kde_points',256, ...
    'compute_kde',true,'rank_tol',128*eps,'support_tol',512*eps, ...
    'integration_abs_tol',2e-7,'integration_rel_tol',2e-5,'kde_max_points',16384);
fn=fieldnames(defaults);
for k=1:numel(fn)
    if ~isfield(opt,fn{k}), opt.(fn{k})=defaults.(fn{k}); end
end
assert(opt.ellipse_points>=16 && opt.quadrature_order>=8 && opt.kde_points>=32 && ...
    opt.rank_tol>0 && opt.rank_tol<1e-3,'RPT:InvalidSeparationOptions', ...
    'Invalid separation resolution or numerical rank tolerance.');
opt.ellipse_points=round(opt.ellipse_points);
opt.quadrature_order=round(opt.quadrature_order);
opt.kde_points=round(opt.kde_points);
end

function [V,lambda,r]=covariance_axes(C,rel_tol)
C=.5*(C+C');
[V,D]=eig(C); [lambda,idx]=sort(real(diag(D))); V=real(V(:,idx));
tol=rel_tol*max(max(lambda),0);
assert(min(lambda)>=-max(tol,64*eps*max(norm(C,'fro'),realmin)), ...
    'RPT:NonPSDCovariance','Covariance is not positive semidefinite.');
lambda(lambda<=tol)=0; r=sum(lambda>0);
end

function ovl=gaussian_overlap_1d(m1,s1,m2,s2)
% Exact normal-density overlap, including point masses. No Statistics Toolbox.
assert(all(isfinite([m1 s1 m2 s2])) && s1>=0 && s2>=0, ...
    'RPT:InvalidNormal','Normal means/scales must be finite; scales nonnegative.');
if s1==0 && s2==0, ovl=double(m1==m2); return; end
if s1==0 || s2==0, ovl=0; return; end
if s1==s2
    ovl=erfc(abs(m2-m1)/(2*sqrt(2)*s1));
else
    ovl=weighted_normal_overlap(0,m1,s1,0,m2,s2);
end
ovl=min(max(real(ovl),0),1);
end

function ov=weighted_normal_overlap(la,m1,s1,lb,m2,s2)
% Integral min(exp(la)*N1,exp(lb)*N2), vectorized over means/log weights.
% Centering/scaling, stable quadratic roots and erfc tails avoid cancellation.
shape=size(la+m1+lb+m2);
la=la+zeros(shape); lb=lb+zeros(shape);
m1=m1+zeros(shape); m2=m2+zeros(shape);
la=la(:); lb=lb(:); d=(m2(:)-m1(:))/max(s1,s2);
a=s1/max(s1,s2); b=s2/max(s1,s2);
logscale=max(la,lb); wa=exp(la-logscale); wb=exp(lb-logscale);
v=zeros(size(d));
if abs(a-b)<=32*eps
    identical=(d==0); v(identical)=min(wa(identical),wb(identical));
    active=~identical;
    t=.5*d(active)+(la(active)-lb(active))./d(active);
    ua=t/a; ub=(t-d(active))/b;
    left_c=(d(active)<0);
    vv=wa(active).*.5.*erfc(ua/sqrt(2)) + ...
       wb(active).*.5.*erfc(-ub/sqrt(2));
    vv(left_c)=wa(active & d<0).*.5.*erfc(-ua(left_c)/sqrt(2)) + ...
       wb(active & d<0).*.5.*erfc(ub(left_c)/sqrt(2));
    v(active)=vv;
else
    A=a^2-b^2;
    B=-2*d*a^2;
    C=d.^2*a^2+2*a^2*b^2*(log(b/a)+la-lb);
    disc=B.^2-4*A*C;
    no_cross=disc<=0;
    v(no_cross)=min(wa(no_cross),wb(no_cross));
    active=~no_cross;
    sg=ones(sum(active),1); sg(B(active)<0)=-1;
    q=-.5*(B(active)+sg.*sqrt(max(disc(active),0)));
    r1=q/A; r2=C(active)./q;
    lo=min(r1,r2); hi=max(r1,r2);
    d0=d(active);
    if A<0
        tails=.5*erfc(-lo/(a*sqrt(2)))+.5*erfc(hi/(a*sqrt(2)));
        mid=normal_interval_mass((lo-d0)/b,(hi-d0)/b);
        v(active)=wa(active).*tails+wb(active).*mid;
    else
        tails=.5*erfc(-(lo-d0)/(b*sqrt(2)))+.5*erfc((hi-d0)/(b*sqrt(2)));
        mid=normal_interval_mass(lo/a,hi/a);
        v(active)=wb(active).*tails+wa(active).*mid;
    end
end
ov=reshape(exp(logscale).*max(v,0),shape);
end

function p=normal_interval_mass(lo,hi)
% Phi(hi)-Phi(lo), avoiding subtraction of two values close to one.
p=zeros(size(lo)); L=hi<=0; R=lo>=0; M=~(L|R);
p(L)=.5*(erfc(-hi(L)/sqrt(2))-erfc(-lo(L)/sqrt(2)));
p(R)=.5*(erfc(lo(R)/sqrt(2))-erfc(hi(R)/sqrt(2)));
p(M)=.5*(erf(hi(M)/sqrt(2))-erf(lo(M)/sqrt(2)));
p=max(p,0);
end

function [ov,err,tail,method]=gaussian_overlap_2d(d,Cc,Ch,Vc,lc,rc,Vh,lh,rh,Vp,lp,opt)
err=0; tail=0;
if rc==0 && rh==0
    ov=double(norm(d)<=opt.support_tol); method='two atoms'; return;
elseif rc~=rh || rc==0 || rh==0
    ov=0; method='mutually singular supports'; return;
elseif rc==1
    vc=Vc(:,lc>0); vh=Vh(:,lh>0);
    same_line=abs(vc(1)*vh(2)-vc(2)*vh(1))<=opt.support_tol && ...
        norm(d'-(vc'*d')*vc)<=opt.support_tol*(1+norm(d));
    if same_line
        ov=gaussian_overlap_1d(0,sqrt(lc(lc>0)),d*vc,sqrt(lh(lh>0)));
        method='common singular line';
    else
        ov=0; method='distinct singular lines';
    end
    return;
end

% Common whitening is invertible and cannot change probability overlap.
W=Vp*diag(1./sqrt(lp));
dw=d*W;
if norm(Cc-Ch,'fro')<=opt.rank_tol*max(1,norm(.5*(Cc+Ch),'fro'))
    ov=erfc(norm(dw)/(2*sqrt(2))); method='equal-covariance analytic'; return;
end
if norm(dw)>opt.support_tol
    ey=dw'/norm(dw); ex=[-ey(2);ey(1)]; R=[ex ey];
else
    R=eye(2);
end
W=W*R; d=d*W;
A=W'*Cc*W; B=W'*Ch*W;
A=.5*(A+A'); B=.5*(B+B');
mx1=-d(1)/2; mx2=d(1)/2; my1=-d(2)/2; my2=d(2)/2;
sx1=sqrt(A(1,1)); sx2=sqrt(B(1,1));
beta1=A(1,2)/A(1,1); beta2=B(1,2)/B(1,1);
% Positive Cholesky pivots are more reliable than subtracting correlated variances.
L1=chol(A,'lower'); L2=chol(B,'lower');
sy1=L1(2,2); sy2=L2(2,2);
fun=@(x) conditional_gaussian_overlap(x,mx1,sx1,my1,beta1,sy1, ...
    mx2,sx2,my2,beta2,sy2);
lo=min(mx1-9*sx1,mx2-9*sx2); hi=max(mx1+9*sx1,mx2+9*sx2);
tail=erfc(9/sqrt(2)); % rigorous bound on omitted overlap outside both 9-sigma ranges
% Two fixed quadrature grids can BOTH miss a thin overlap band. Split adaptive
% integration at both marginal scales and at the conditional-mean crossing.
% quadrature_order controls the starting split density; quadgk refines further.
levels=sort(unique([-9 -4 -2 -1 0 1 2 4 9, ...
    linspace(-4,4,max(3,round(opt.quadrature_order/6)))]));
waypoints=[mx1+levels*sx1, mx2+levels*sx2];
db=beta1-beta2;
if db~=0
    cross_x=(my2-beta2*mx2-my1+beta1*mx1)/db;
    cross_width=(sy1+sy2)/abs(db);
    waypoints=[waypoints, cross_x+levels*cross_width];
end
waypoints=sort(unique(waypoints(isfinite(waypoints) & waypoints>lo & waypoints<hi)));
[ov,err0]=quadgk(fun,lo,hi,'Waypoints',waypoints, ...
    'AbsTol',opt.integration_abs_tol,'RelTol',opt.integration_rel_tol, ...
    'MaxIntervalCount',max(650,32*opt.quadrature_order));
err=err0+tail; method='conditional analytic + split adaptive Gauss-Kronrod';
ov=min(max(real(ov),0),1);
end

function y=conditional_gaussian_overlap(x,mx1,sx1,my1,beta1,sy1,mx2,sx2,my2,beta2,sy2)
la=-log(sqrt(2*pi)*sx1)-.5*((x-mx1)/sx1).^2;
lb=-log(sqrt(2*pi)*sx2)-.5*((x-mx2)/sx2).^2;
m1=my1+beta1*(x-mx1); m2=my2+beta2*(x-mx2);
y=weighted_normal_overlap(la,m1,sy1,lb,m2,sy2);
end

function S=ellipse_jaccard_separation(mc,Cc,mh,Ch,n)
% Polygon approximation to ellipse intersection/union, not probability overlap.
if isequal(mc,mh) && norm(Cc-Ch,'fro')<=128*eps*max(1,norm(Cc,'fro'))
    S=0; return;
end
p=ellipse_polygon(mc,Cc,n); q=ellipse_polygon(mh,Ch,n);
if any(max(p,[],1)<min(q,[],1)) || any(max(q,[],1)<min(p,[],1))
    S=1; return;
end
ap=convex_polygon_area(p); aq=convex_polygon_area(q);
if ap==0 || aq==0, S=NaN; return; end
% This triangle/spectral-norm bound is sufficient for full-ellipse containment;
% it never declares containment merely because sampled vertices happen to fit.
if ellipse_contains(mc,Cc,mh,Ch), S=1-aq/ap; return; end
if ellipse_contains(mh,Ch,mc,Cc), S=1-ap/aq; return; end
r=p;
for k=1:size(q,1)
    a=q(k,:); b=q(mod(k,size(q,1))+1,:); edge=b-a;
    if isempty(r), break; end
    prev=r([end 1:end-1],:);
    dc=edge(1)*(r(:,2)-a(2))-edge(2)*(r(:,1)-a(1));
    dp=edge(1)*(prev(:,2)-a(2))-edge(2)*(prev(:,1)-a(1));
    tol=64*eps*norm(edge)*max(1,max(abs([r(:);a(:);b(:)])));
    inside=dc>=-tol; wasinside=dp>=-tol; crossing=xor(inside,wasinside);
    count=double(inside)+double(crossing); finish=cumsum(count);
    out=zeros(sum(count),2);
    out(finish(inside),:)=r(inside,:);
    t=dp(crossing)./(dp(crossing)-dc(crossing));
    out(finish(crossing)-double(inside(crossing)),:)=prev(crossing,:)+ ...
        bsxfun(@times,t,r(crossing,:)-prev(crossing,:));
    r=out;
end
ai=min([convex_polygon_area(r),ap,aq]);
S=min(max(1-ai/(ap+aq-ai),0),1);
end

function yes=ellipse_contains(mc,Cc,mh,Ch)
[V,D]=eig(.5*(Cc+Cc')); [U,E]=eig(.5*(Ch+Ch'));
scale=sqrt(max(diag(D),realmin));
offset=bsxfun(@rdivide,(mh-mc)*V,scale'*sqrt(-2*log(.05)));
B=bsxfun(@rdivide,V'*U*diag(sqrt(max(diag(E),0))),scale);
yes=norm(offset)+norm(B,2)<=1-128*eps;
end

function p=ellipse_polygon(mu,C,n)
[V,D]=eig(.5*(C+C')); t=(0:n-1)'*(2*pi/n);
p=bsxfun(@plus,[cos(t) sin(t)]*diag(sqrt(max(diag(D),0)*(-2*log(.05))))*V',mu);
if signed_polygon_area(p)<0, p=flipud(p); end
end

function a=convex_polygon_area(p)
if size(p,1)<3, a=0; else, a=abs(signed_polygon_area(p)); end
end

function a=signed_polygon_area(p)
p=bsxfun(@minus,p,mean(p,1)); q=p([2:end 1],:);
a=.5*sum(p(:,1).*q(:,2)-q(:,1).*p(:,2));
end

function [ov,grid,fc,fh,info]=projected_kde_overlap(c,h,sc,sh,opt)
% Gaussian KDE with Silverman's bandwidth and linear binning + FFT convolution.
% Smoothing estimates a density; it is not added physical observation noise.
grid=[]; fc=[]; fh=[];
info=struct('points',0,'step',NaN,'bandwidths',[0 0], ...
    'tail_bound',0,'resolved',true);
if sc==0 && sh==0, ov=double(c(1)==h(1)); return; end
if sc==0 || sh==0, ov=0; return; end
bc=1.06*sc*numel(c)^(-.2); bh=1.06*sh*numel(h)^(-.2);
info.bandwidths=[bc bh];
if isequal(c,h), ov=1; return; end
if max(c)+8*bc<min(h)-8*bh || max(h)+8*bh<min(c)-8*bc
    ov=0; info.tail_bound=erfc(8/sqrt(2)); return;
end
lo=min(min(c)-8*bc,min(h)-8*bh); hi=max(max(c)+8*bc,max(h)+8*bh);
needed=ceil(2*(hi-lo)/min(bc,bh))+1;
n=2^nextpow2(max(opt.kde_points,needed));
if n>opt.kde_max_points
    % Refuse an unresolved narrow KDE rather than silently creating a density floor.
    ov=NaN; info.points=n; info.resolved=false; return;
end
grid=linspace(lo,hi,n)'; dx=grid(2)-grid(1);
fc=binned_gaussian_kde(c,grid,bc); fh=binned_gaussian_kde(h,grid,bh);
ov=min(max(trapz(grid,min(fc,fh)),0),1);
info.points=n; info.step=dx; info.tail_bound=2*erfc(8/sqrt(2));
end

function f=binned_gaussian_kde(z,grid,bw)
n=numel(grid); dx=grid(2)-grid(1);
u=(z(:)-grid(1))/dx+1; k=floor(u); a=u-k;
mass=accumarray([k;k+1],[1-a;a],[n 1])/numel(z);
radius=ceil(8*bw/dx); offsets=(-radius:radius)'*dx;
kernel=exp(-.5*(offsets/bw).^2)/(bw*sqrt(2*pi));
nf=2^nextpow2(n+numel(kernel)-1);
allf=real(ifft(fft(mass,nf).*fft(kernel,nf)));
f=max(allf(radius+(1:n)),0);
f=f/trapz(grid,f);
end


function run_test3_self_checks()
% Deterministic checks of the statistical definitions, independent of CoolProp.
% A lost null-space direction, a fixed dimensional ridge, an ellipse/S confusion,
% or a missed Gaussian tail makes at least one check below fail.
opt = struct('compute_kde', true, 'kde_points', 256, ...
    'quadrature_order', 24, 'ellipse_points', 128);
Q = sqrt(3/2)*[-1 0; 1 0; 0 -1; 0 1]; % sample mean=0, covariance=I
a = ensemble_separability(Q, Q+[2 0], opt);
expected = 1-erfc(1/sqrt(2));
assert(abs(a.S_lda-expected)<2e-12 && ...
       abs(a.S_gaussian2d-expected)<2e-12, ...
       'RPT:CheckEqualGaussian', 'Equal-covariance analytic overlap failed.');
assert(abs(a.Pe_lda-.5*(1-a.S_lda))<eps && ...
       abs(a.Pe_gaussian2d-.5*(1-a.S_gaussian2d))<eps);
assert(abs(gaussian_overlap_1d(0,1,0,2)-.6773254311652315)<2e-12);
assert(abs(gaussian_overlap_1d(0,1,2,1)-erfc(1/sqrt(2)))<2e-12);
assert(gaussian_overlap_1d(1,0,1,0)==1 && ...
       gaussian_overlap_1d(1,0,2,0)==0 && ...
       gaussian_overlap_1d(1,0,1,1)==0);

b = ensemble_separability(Q+[2 0],Q,opt);
names = {'S_ellipse','S_lda','S_gaussian2d','S_kde'};
for j=1:numel(names)
    assert(abs(a.(names{j})-b.(names{j}))<2e-7, ...
        'RPT:CheckLabelSymmetry','Class-label symmetry failed.');
end
D = diag([1e6 1e-3]);
b = ensemble_separability(Q*D+[123 -7],(Q+[2 0])*D+[123 -7],opt);
for j=1:numel(names)
    assert(abs(a.(names{j})-b.(names{j}))<2e-7, ...
        'RPT:CheckUnits','Changing feature units changed separability.');
end

a = ensemble_separability(Q,Q,opt);
assert(abs(a.S_ellipse)<2e-10 && a.S_lda==0 && ...
    a.S_gaussian2d==0 && a.S_kde==0, ...
    'RPT:CheckIdentical','Identical distributions must have zero separation.');
a = ensemble_separability(Q,Q+[20 0],opt);
assert(a.S_ellipse==1,'RPT:CheckDisjointEllipses','Disjoint ellipses failed.');
a = ensemble_separability(Q,2*Q,opt);
assert(abs(a.S_ellipse-.75)<2e-10, ...
    'RPT:CheckContainedEllipses','Ellipse area ratio must survive containment.');

a = ensemble_separability(Q,Q*diag([2 .5]),opt);
assert(isnan(a.S_lda) && isnan(a.S_kde) && ~a.diagnostics.lda_defined, ...
    'RPT:CheckSameMean','Equal means do not define a Fisher LDA direction.');
assert(a.S_gaussian2d>.1 && a.S_gaussian2d<1 && ...
    a.diagnostics.gaussian2d_error<1e-5, ...
    'RPT:CheckCovarianceContrast','2D Gaussian covariance contrast was lost.');
for narrow=[.01 .001 .00001]
    a=ensemble_separability(Q*diag([narrow 1]),Q,opt);
    expected=1-gaussian_overlap_1d(0,narrow,0,1);
    assert(abs(a.S_gaussian2d-expected)<2e-7, ...
        'RPT:CheckNarrowGaussian','Both quadrature grids missed narrow probability mass.');
end
R45=[1 -1;1 1]/sqrt(2);
for narrow=[.01 .001 .00001]
    a=ensemble_separability(Q*diag([narrow 1])*R45', ...
        Q*diag([1 narrow])*R45',opt);
    assert(abs(a.S_gaussian2d-(1-4*atan(narrow)/pi))<2e-7, ...
        'RPT:CheckCrossedBands','Conditional integration missed thin crossed bands.');
end

% Identical Gaussian moments do not define a direction for non-Gaussian KDEs.
a=ensemble_separability(Q,Q*R45,opt);
assert(a.S_gaussian2d<2e-12 && a.S_lda==0 && isnan(a.S_kde), ...
    'RPT:CheckUndefinedKDE','An arbitrary direction must not masquerade as Fisher KDE.');

t = (-5:5)'; L = [t zeros(size(t))];
a = ensemble_separability(L,L+[0 .01],opt);
assert(a.rank_c==1 && a.rank_h==1 && isnan(a.S_ellipse) && ...
    a.S_lda==1 && a.S_gaussian2d==1 && a.S_kde==1, ...
    'RPT:CheckParallelLines','Distinct singular lines must remain separable.');
a = ensemble_separability(L,[zeros(size(t)) t],opt);
assert(a.S_gaussian2d==1, ...
    'RPT:CheckIntersectingLines','Distinct singular supports share zero mass.');
a = ensemble_separability(L,L+[2 0],opt);
assert(abs(a.S_gaussian2d-(1-erfc(1/(sqrt(2)*std(t)))))<2e-12);

C = repmat([10 2],12,1);
a = ensemble_separability(C,C,opt);
assert(a.S_lda==0 && a.S_gaussian2d==0 && a.S_kde==0 && isnan(a.S_ellipse));
a = ensemble_separability(C,C+[.001 0],opt);
assert(a.S_lda==1 && a.S_gaussian2d==1 && a.S_kde==1 && isnan(a.S_ellipse));
a = ensemble_separability(C,[t+10 zeros(size(t))+2],opt);
assert(a.S_gaussian2d==1,'RPT:CheckAtom','An atom and a continuous law share zero mass.');

R = [cos(.4) -sin(.4); sin(.4) cos(.4)];
C = Q*diag([1 1e-5])*R'; H=C+([0 5e-6]*R');
a = ensemble_separability(C,H,opt);
b = ensemble_separability(C*diag([1e5 1e-3]),H*diag([1e5 1e-3]),opt);
assert(a.rank_c==2 && a.rank_h==2 && ...
    abs(a.S_lda-b.S_lda)<2e-5 && abs(a.S_gaussian2d-b.S_gaussian2d)<2e-5, ...
    'RPT:CheckThinCloud','A thin but resolved covariance was discarded.');

H = Q*[1.3 .2; .1 .7]+[1.2 .3];
a = ensemble_separability(Q,H,opt);
assert(a.S_gaussian2d+2e-6>=a.S_lda && a.diagnostics.gaussian2d_error<1e-5, ...
    'RPT:CheckProjectionInformation','Projection cannot improve Gaussian Bayes separation.');
opt.kde_points=512;
b = ensemble_separability(Q,H,opt);
assert(abs(a.S_kde-b.S_kde)<3e-3, ...
    'RPT:CheckKDEResolution','KDE overlap did not converge with grid refinement.');

% 99% separability: this overlap must survive numerical tail handling.
a = ensemble_separability(Q,Q+[5.151658607097801 0],opt);
assert(abs(a.S_gaussian2d-.99)<2e-12 && abs(a.S_lda-.99)<2e-12, ...
    'RPT:CheckRareOverlap','Rare Gaussian overlap was rounded away.');
fprintf('test3 statistical self-checks passed.\n');
end


function fig = render_separability_comparison(map, lang, cfg)
% Display one four-method comparison; this function never writes a file.
% Each S_* field is ny-by-nx in [0,1], with NaN for undefined cells.
% x/y are the actual parameter coordinates, already in display units.
if nargin < 2 || isempty(lang), lang = 'EN'; end
if nargin < 3 || isempty(cfg), cfg = struct(); end
lang = upper(char(lang));
fields = {'S_ellipse', 'S_lda', 'S_gaussian2d', 'S_kde'};
required = [{'name', 'kind', 'x', 'y', 'xlabel_cn', 'xlabel_en', ...
    'ylabel_cn', 'ylabel_en'}, fields];
for k = 1:numel(required)
    assert(isfield(map, required{k}), 'RPT:InvalidSeparabilityMap', ...
        'The comparison map is missing field %s.', required{k});
end
assert(isnumeric(map.x) && isreal(map.x) && isvector(map.x) && ...
    ~isempty(map.x) && all(isfinite(map.x(:))) && ...
    isnumeric(map.y) && isreal(map.y) && isvector(map.y) && ...
    ~isempty(map.y) && all(isfinite(map.y(:))), ...
    'RPT:InvalidSeparabilityMap', 'Map coordinates must be finite real vectors.');
[x, ix] = sort(double(map.x(:)'));
[y, iy] = sort(double(map.y(:)));
assert(all(diff(x) > 0) && all(diff(y) > 0), ...
    'RPT:InvalidSeparabilityMap', 'Map coordinates must be distinct.');
[X, Y] = meshgrid(x, y);
nx = numel(x); ny = numel(y);
for k = 1:numel(fields)
    values = map.(fields{k});
    assert(isnumeric(values) && isreal(values) && ...
        isequal(size(values), [ny nx]), ...
        'RPT:InvalidSeparabilityMap', '%s must be ny-by-nx.', fields{k});
    finite_values = values(~isnan(values));
    assert(all(isfinite(finite_values)) && ...
        all(finite_values >= -1e-12 & finite_values <= 1+1e-12), ...
        'RPT:InvalidSeparabilityMap', '%s must be in [0,1] or NaN.', fields{k});
end

position = [60 60 1200 850];
if isfield(cfg, 'figure_position') && ~isempty(cfg.figure_position)
    assert(isnumeric(cfg.figure_position) && numel(cfg.figure_position) == 4 && ...
        all(isfinite(cfg.figure_position(:))) && ...
        all(cfg.figure_position(3:4) > 0), ...
        'RPT:InvalidFigurePosition', 'figure_position must be [left bottom width height].');
    position = double(reshape(cfg.figure_position, 1, 4));
end
fig = figure('Name', sprintf('%s | %s', char(map.name), lang), ...
    'NumberTitle', 'off', 'Visible', 'on', 'Color', 'white', ...
    'Units', 'pixels', 'Position', position);
layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'loose');

cn_titles = {'(a) 95% 协方差椭圆：面积启发式评分', ...
    '(b) LDA 投影：一维高斯重叠', ...
    '(c) 二维高斯重叠', '(d) LDA 投影：一维 KDE 重叠'};
en_titles = {'(a) 95% covariance ellipses: area heuristic', ...
    '(b) LDA projection: 1D Gaussian overlap', ...
    '(c) 2D Gaussian overlap', '(d) LDA projection: 1D KDE overlap'};
gray = [.82 .83 .85];
levels = [50 90 99];
contour_colors = {[1 1 1], [.10 .10 .10], [.65 .05 .35]};
contour_styles = {'-', '--', ':'};
is_reference_map = strncmpi(char(map.kind), 'B_', 2);
show_reference = is_reference_map && isfield(map, 'reference_xy') && ...
    isnumeric(map.reference_xy) && numel(map.reference_xy) == 2 && ...
    all(isfinite(map.reference_xy(:)));

for k = 1:4
    ax = nexttile(layout, k);
    hold(ax, 'on');
    set(ax, 'Color', gray, 'YDir', 'normal', 'FontName', 'Times New Roman', ...
        'FontSize', 10.5, 'FontWeight', 'bold', 'LineWidth', 1.05, ...
        'TickDir', 'out', 'Box', 'on', 'Layer', 'top', 'SortMethod', 'childorder', ...
        'XColor', [.18 .22 .26], 'YColor', [.18 .22 .26]);
    raw = double(map.(fields{k}));
    raw = raw(iy, ix);
    Z = 100*min(1, max(0, raw));
    Z(isnan(raw)) = NaN;
    sep_draw_flat_cells(ax, x, y, Z);
    colormap(ax, parula(256));
    caxis(ax, [0 100]);
    xlim(ax, sep_plot_limits(x));
    ylim(ax, sep_plot_limits(y));
    [txt, intrp, fnt] = figure_label(choose_label(lang, ...
        map.xlabel_cn, map.xlabel_en), lang);
    xlabel(ax, txt, 'Interpreter', intrp, 'FontName', fnt, ...
        'FontSize', 12, 'FontWeight', 'bold');
    [txt, intrp, fnt] = figure_label(choose_label(lang, ...
        map.ylabel_cn, map.ylabel_en), lang);
    ylabel(ax, txt, 'Interpreter', intrp, 'FontName', fnt, ...
        'FontSize', 12, 'FontWeight', 'bold');
    [txt, intrp, fnt] = figure_label(choose_label(lang, ...
        cn_titles{k}, en_titles{k}), lang);
    title(ax, txt, 'Interpreter', intrp, 'FontName', fnt, ...
        'FontSize', 12.5, 'FontWeight', 'bold');

    cb = colorbar(ax);
    set(cb, 'Ticks', 0:20:100, 'FontName', 'Times New Roman', ...
        'FontSize', 10, 'FontWeight', 'bold', 'LineWidth', 1);
    if k == 1
        cb_cn = '95% 椭圆面积评分 (%)';
        cb_en = '95% ellipse-area score (%)';
    else
        cb_cn = '可分离度 S = 1 - 2P_e (%)';
        cb_en = 'Separability S = 1 - 2P_e (%)';
    end
    [txt, intrp, fnt] = figure_label(choose_label(lang, cb_cn, cb_en), lang);
    set(cb.Label, 'String', txt, 'Interpreter', intrp, 'FontName', fnt, ...
        'FontSize', 11, 'FontWeight', 'bold');

    % Every legend entry corresponds to a contour actually drawn. No
    % extrapolated contours are manufactured when a threshold is absent.
    key_handles = gobjects(0, 1);
    key_labels = cell(0, 1);
    finite_z = Z(isfinite(Z));
    if nx >= 2 && ny >= 2 && ~isempty(finite_z)
        for j = 1:numel(levels)
            level = levels(j);
            if level <= min(finite_z) || level >= max(finite_z), continue; end
            [cc, hh] = contour(ax, X, Y, Z, [level level], ...
                'LineColor', contour_colors{j}, 'LineStyle', contour_styles{j}, ...
                'LineWidth', 1.35);
            if isempty(cc)
                set(hh, 'HandleVisibility', 'off');
                continue;
            end
            clabel(cc, hh, 'FontName', 'Times New Roman', 'FontSize', 8.5, ...
                'FontWeight', 'bold', 'Color', contour_colors{j}, 'LabelSpacing', 260);
            key_handles(end+1, 1) = hh; %#ok<AGROW>
            if k == 1
                key_labels{end+1, 1} = choose_label(lang, ...
                    sprintf('面积评分 %g%%', level), sprintf('Area score %g%%', level)); %#ok<AGROW>
            else
                key_labels{end+1, 1} = sprintf('S = %g%%', level); %#ok<AGROW>
            end
        end
    end
    if isempty(finite_z)
        [txt, intrp, fnt] = figure_label(choose_label(lang, ...
            '此网格无有效数值', 'No defined values on this grid'), lang);
        text(ax, .5, .5, txt, 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'Interpreter', intrp, 'FontName', fnt, 'FontSize', 12, ...
            'Color', [.23 .25 .28], 'FontWeight', 'bold');
    end
    if show_reference
        ref_cn = '固定参考点'; ref_en = 'Fixed reference point';
        if isfield(map, 'reference_label_cn') && ~isempty(map.reference_label_cn)
            ref_cn = map.reference_label_cn;
        end
        if isfield(map, 'reference_label_en') && ~isempty(map.reference_label_en)
            ref_en = map.reference_label_en;
        end
        hh = plot(ax, map.reference_xy(1), map.reference_xy(2), 'p', ...
            'MarkerSize', 12, 'MarkerFaceColor', 'white', ...
            'MarkerEdgeColor', [.08 .08 .08], 'LineWidth', 1.5, 'Clipping', 'off');
        key_handles(end+1, 1) = hh; %#ok<AGROW>
        key_labels{end+1, 1} = choose_label(lang, ref_cn, ref_en); %#ok<AGROW>
    end
    if ~isempty(key_handles)
        [txt, intrp, fnt] = figure_label(key_labels, lang);
        legend(ax, key_handles, txt, 'Location', 'southoutside', ...
            'Orientation', 'horizontal', 'NumColumns', 2, ...
            'Box', 'off', 'Interpreter', intrp, 'FontName', fnt, ...
            'FontSize', 9.5, 'AutoUpdate', 'off');
    end
end

heading = sep_comparison_heading(map, lang);
notes = cell(0, 1);
note_field = choose_label(lang, 'note_cn', 'note_en');
if isfield(map, note_field) && ~isempty(map.(note_field))
    notes = sep_note_lines(map.(note_field));
end
semantics = choose_label(lang, ...
    '灰色：未定义。椭圆评分为面积指标；其余三图等先验 S = 1 - 2P_e，S ≥ 99% 表示近乎完全可分。', ...
    'Gray: undefined. Ellipse area is heuristic; other panels: equal-prior S = 1 - 2P_e; S >= 99% means near-complete separation.');
[txt, intrp, fnt] = figure_label([{heading}; notes; {semantics}], lang);
title(layout, txt, 'Interpreter', intrp, 'FontName', fnt, ...
    'FontSize', 12.5, 'FontWeight', 'bold');
drawnow;
end

function sep_draw_flat_cells(ax, x, y, Z)
% One flat-colored cell per evaluated coordinate; no interpolated shading.
% Explicit cell edges preserve the spacing of nonuniform parameter grids.
xe = sep_cell_edges(x); ye = sep_cell_edges(y);
[Xe, Ye] = meshgrid(xe, ye);
Zpad = [Z Z(:, end); Z(end, :) Z(end, end)];
alpha = double(isfinite(Zpad));
Zpad(~isfinite(Zpad)) = 0;
surface(ax, Xe, Ye, zeros(size(Xe)), Zpad, ...
    'FaceColor', 'flat', 'EdgeColor', 'none', 'CDataMapping', 'scaled', ...
    'FaceAlpha', 'flat', 'AlphaData', alpha, 'AlphaDataMapping', 'none', ...
    'HandleVisibility', 'off');
view(ax, 2);
end

function edges = sep_cell_edges(v)
v = v(:)';
if numel(v) == 1
    edges = sep_plot_limits(v);
else
    % Boundary cells end at the evaluated parameter domain, not beyond it.
    edges = [v(1), v(1:end-1) + diff(v)/2, v(end)];
end
end

function limits = sep_plot_limits(v)
v = v(:);
if numel(v) == 1
    padding = max(.5, .02*abs(v(1)));
    limits = [v(1)-padding v(1)+padding];
else
    limits = [v(1) v(end)];
end
end

function heading = sep_comparison_heading(map, lang)
switch lower(char(map.kind))
    case 'a_phi'
        cn = '图 3 | 方案 A：孔隙度与气体饱和度可分离度';
        en = 'Figure 3 | Scheme A: porosity and gas-saturation separability';
    case 'a_tp'
        cn = '图 4 | 方案 A：温度与压力可分离度';
        en = 'Figure 4 | Scheme A: temperature and pressure separability';
    case 'b_phi'
        cn = '图 5 | 方案 B：孔隙度与气体饱和度可分离度';
        en = 'Figure 5 | Scheme B: porosity and gas-saturation separability';
    case 'b_tp'
        cn = '图 6 | 方案 B：温度与压力可分离度';
        en = 'Figure 6 | Scheme B: temperature and pressure separability';
    otherwise
        cn = char(map.name); en = char(map.name);
end
heading = choose_label(lang, cn, en);
end

function lines = sep_note_lines(note)
if ischar(note)
    lines = regexp(note, '\r\n|\n|\r', 'split')';
elseif isstring(note)
    lines = cellstr(note(:));
elseif iscell(note)
    lines = note(:);
else
    error('RPT:InvalidMapNote', 'Comparison notes must be text or a cell array of text.');
end
lines = lines(~cellfun(@isempty, lines));
end

function E = covariance_ellipse_points(X, level, n)
% Return an n-by-2 covariance contour in X's original physical units.
% The Gaussian 2D probability contour has radius sqrt(-2*log(1-level)).
% Rank-one and rank-zero clouds produce a line and a repeated point.
% No variance floor, chi2inv, or other Statistics Toolbox call is used.
if nargin < 2 || isempty(level), level = .95; end
if nargin < 3 || isempty(n), n = 361; end
assert(isnumeric(X) && isreal(X) && ismatrix(X) && size(X, 2) == 2, ...
    'RPT:InvalidEllipseSamples', 'Ellipse samples must be a real N-by-2 matrix.');
assert(isnumeric(level) && isscalar(level) && isfinite(level) && level > 0 && level < 1, ...
    'RPT:InvalidEllipseLevel', 'The ellipse level must lie strictly between zero and one.');
assert(isnumeric(n) && isscalar(n) && isfinite(n) && n >= 3 && n == fix(n), ...
    'RPT:InvalidEllipsePointCount', 'The ellipse needs an integer number of points n >= 3.');
X = double(X(all(isfinite(X), 2), :));
if isempty(X)
    E = NaN(n, 2);
    return;
end
mu = mean(X, 1);
if size(X, 1) == 1
    E = repmat(mu, n, 1);
    return;
end
centered = bsxfun(@minus, X, mu);
scale = max(abs(centered), [], 1);
scale(scale == 0) = 1;
% Independently scale coordinates before eigendecomposition so changing
% impedance units cannot erase a real, small variance on the ratio axis.
U = bsxfun(@rdivide, centered, scale);
C = (U'*U)/(size(U, 1)-1);
C = (C+C')/2;
[V, D] = eig(C);
lambda = max(real(diag(D)), 0); % Clip negative roundoff only.
% n distinct periodic points also gives a valid three-vertex contour at n=3.
theta = (0:n-1)'*(2*pi/n);
radius = sqrt(-2*log1p(-level));
unit_circle = [cos(theta), sin(theta)];
curve_scaled = radius*unit_circle*diag(sqrt(lambda))*V';
E = bsxfun(@plus, bsxfun(@times, curve_scaled, scale), mu);
end

function h = add_covariance_ellipse(ax, X, color, level)
% Add a dashed covariance contour and a mean marker to an existing axes.
% Keep existing traces and legend membership; do not open another figure.
if nargin < 4 || isempty(level), level = .95; end
E = covariance_ellipse_points(X, level, 361);
h = gobjects(0, 1);
if isempty(E) || ~any(all(isfinite(E), 2)), return; end
Xvalid = double(X(all(isfinite(X), 2), :));
mu = mean(Xvalid, 1);
was_held = ishold(ax);
hold(ax, 'on');
h(1, 1) = plot(ax, [E(:, 1); E(1, 1)], [E(:, 2); E(1, 2)], '--', 'Color', color, ...
    'LineWidth', 1.7, 'HandleVisibility', 'off');
h(2, 1) = plot(ax, mu(1), mu(2), '+', 'Color', color, ...
    'MarkerSize', 8, 'LineWidth', 1.8, 'HandleVisibility', 'off');
if ~was_held, hold(ax, 'off'); end
end
