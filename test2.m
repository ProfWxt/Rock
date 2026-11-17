% MATLAB ROCK PHYSICS - STANDALONE BILINGUAL PANELS (2026-11-11)
% Figure 1: (a,b,d,e,f,g) independent 2.90 x 3.60 inch PORTRAIT panels.
% Figure 1(c) and Figure 2(a-d): independent 6.40 x 3.60 inch LANDSCAPE.
% Every independent panel has the same height: 2/5 of a 16 x 9-inch PPT slide.
% Figure 1 and Figure 2 use separate output folders; CN and EN files coexist.
% Keep the letter identifiers (a), (b), ...; NO panel titles or metric tags.
% Figure 1(c) has NO zoom inset; all main axes have 11.5-pt PDF labels.
% TeX-label scaling is calibrated against Fig1_a_EN.pdf supplied 2026-11-11.
% Figure 1(a,b,d,e,f,g) have identical plot widths/heights, not only pages.
% Physical conditions, axis labels, numerical contours, and essential keys stay.
% Scientific calculations, gas properties, and sampling are unchanged.
% MATLAB + CoolProp required; SVG/PDF vector export uses painters.
clear; close all; clc;
% Scientific text specifies its interpreter at each call site; do not
% change MATLAB's session-wide default interpreter or font preferences.
check_figure_fonts();
% 部分饱和采用 Gassmann-Hill 斑块饱和；气体物性必须由 CoolProp 提供。
rng(2026, 'twister'); % 固定随机种子以保证蒙特卡洛结果完全可复现

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
phi_fixed  = 0.10;   % 固定孔隙度 10%
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
T_rpt = 60;         % 温度基准 [degC]
P_rpt = 35e6;       % 孔隙压力基准 [Pa]

N_phi = 100; N_sg = 100;
phi_vec = linspace(0.05, 0.30, N_phi);
Sg_vec  = linspace(0, 1, N_sg);
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
gas_lut = build_coolprop_lut(linspace(35, 125, 45), ...
    linspace(0.75e6, 46e6, 55));

% P-T 热图：保持 90x90 分辨率和独立的线性色标。
P_pt = linspace(1e6, 35e6, 90);
T_pt = linspace(40, 120, 90);
[PP_pt, TT_pt] = meshgrid(P_pt, T_pt);
phi_pt_grid = phi_fixed * ones(size(PP_pt));
sg_pt_grid  = 0.50 * ones(size(PP_pt));
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
%  All changes below affect layout/typography ONLY, not the forward model.
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
% Evaluate RPT mesh and sensitivity curves only once for both languages.
A.V = precompute_visual_curves(A);

%% ========================================================================
%  Figure 2: paired parameter ensembles, not measurement noise.
%  This random sampling is executed ONCE so the two language versions match.
%  ========================================================================
N_mc = 3500;
phi_dense_rpt = linspace(0.05, 0.30, 60);
[Ip_sw_line, VpVs_sw_line] = calc_rpt_lut( ...
    interp1(phi_dem_axis, K_dem, phi_dense_rpt, 'pchip'), ...
    interp1(phi_dem_axis, G_dem, phi_dense_rpt, 'pchip'), ...
    phi_dense_rpt, zeros(size(phi_dense_rpt)), T_rpt, P_rpt, ...
    'Methane', salinity, K_solid, rho_solid, gas_lut);

% Row 1: vary ONLY porosity and gas saturation. T and P are STRICTLY fixed.
phi_mc1 = min(max(0.15 + 0.025*randn(N_mc, 1), 0.05), 0.30);
Sg_mc1  = min(max(0.60 + 0.120*randn(N_mc, 1), 0.05), 1.00);
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
Sg_mc2  = ones(N_mc, 1);
T_mc2   = 40 + 80*rand(N_mc, 1);
P_mc2   = (1 + 34*rand(N_mc, 1))*1e6;
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

%% ---------- Separate panels: two folders, PDF + SVG for each panel ----------
% Never generate a 2x4 or 2x2 composite canvas. Every output is one figure.
script_path = mfilename('fullpath');
if isempty(script_path)
    script_dir = pwd;
else
    script_dir = fileparts(script_path);
end
out_root = fullfile(script_dir, 'rockphysics_individual_panels');
fig1_dir = fullfile(out_root, 'Fig1_ElasticContrasts');
fig2_dir = fullfile(out_root, 'Fig2_RPT_LDA');
if ~exist(fig1_dir, 'dir'), mkdir(fig1_dir); end
if ~exist(fig2_dir, 'dir'), mkdir(fig2_dir); end

% 16 x 9 inch PowerPoint reference: 2/5 of slide height = 3.60 inches.
% Portrait: 2.90 x 3.60; landscape: 6.40 x 3.60 (16:9).
% Both page dimensions are expressed in INCHES, not normalized pixels.
slide_height_in = 9;
panel_height_in = (2/5)*slide_height_in;
fig1_portrait_size_in = [2.90, panel_height_in];
landscape_size_in = [(16/9)*panel_height_in, panel_height_in];
assert(fig1_portrait_size_in(1) < panel_height_in && ...
       landscape_size_in(1) > panel_height_in, ...
    'RPT:PanelAspectRatio', ...
    'Figure 1(a,b,d,e,f,g) must be portrait; Figure 1(c) and Figure 2 must be landscape.');
lang_list = {'CN', 'EN'};
fig1_panels = 'abcdefg';
fig2_panels = 'abcd';
for il = 1:numel(lang_list)
    lang = lang_list{il};
    for k = 1:numel(fig1_panels)
        id = fig1_panels(k);
        fig = render_ideal_panel(A, lang, id, ...
            fig1_portrait_size_in, landscape_size_in);
        export_bilingual_figure(fig, fig1_dir, sprintf('Fig1_%s_%s', id, lang));
        close(fig);
    end
    for k = 1:numel(fig2_panels)
        id = fig2_panels(k);
        fig = render_mc_panel(M, lang, k, id, landscape_size_in);
        export_bilingual_figure(fig, fig2_dir, sprintf('Fig2_%s_%s', id, lang));
        close(fig);
    end
end
fprintf(['Saved separate PDF/SVG panels:\n' ...
    '  Figure 1(a,b,d,e,f,g) = %.2f x %.2f inches (portrait)\n' ...
    '  Figure 1(c) and Figure 2(a-d) = %.2f x %.2f inches (landscape)\n' ...
    '  %s\n  %s\n'], ...
    fig1_portrait_size_in(1), fig1_portrait_size_in(2), ...
    landscape_size_in(1), landscape_size_in(2), fig1_dir, fig2_dir);

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
fig = figure('Name', fig_name, 'NumberTitle', 'off', 'Visible', 'off', ...
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
        add_portrait_condition_header(fig, 'T = 60 °C, P = 35 MPa');
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
        xlim(ax, [1 35]); ylim(ax, [40 120]);
        style_colorbar(ax, lang, 'RMS 相对差异 (%)', ...
            'RMS relative contrast (%)');
        add_portrait_condition_header(fig, 'φ = 10%, S_g = 50%');
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
    xlim(ax, [1 35]);
    add_small_condition_note(ax, {'φ = 10%', 'T = 60 °C'}, 'southeast');
else
    apply_axis_labels(ax, lang, ...
        '温度 T (°C)', 'Temperature T (°C)', ...
        'RMS 相对差异 (%)', 'RMS relative contrast (%)');
    xlim(ax, [40 120]);
    add_small_condition_note(ax, {'φ = 10%', 'P = 35 MPa'}, 'southeast');
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
        xticks(ax, [1 10 20 30]);
        yticks(ax, [40 60 80 100 120]);
    else
        xticks(ax, [5 10 15 20 25 30]);
        yticks(ax, [0 20 40 60 80 100]);
    end
elseif id == 'f'
    xticks(ax, [1 10 20 30 35]);
elseif id == 'g'
    xticks(ax, [40 60 80 100 120]);
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
    cn_note = '固定 T = 60 °C、P = 35 MPa；仅 φ、S_g 变化';
    en_note = 'Fixed T = 60 °C, P = 35 MPa; variable φ and S_g';
else
    Ip_c = M.Ip_c_r2; R_c = M.VpVs_c_r2;
    Ip_h = M.Ip_h_r2; R_h = M.VpVs_h_r2;
    brine = false;
    cn_note = '固定 φ = 10%、S_g = 100%；仅 T、P 变化';
    en_note = 'Fixed φ = 10%, S_g = 100%; variable T and P';
end
Ip_c = Ip_c(:); R_c = R_c(:);
Ip_h = Ip_h(:); R_h = R_h(:);
if mod(panel_index, 2) == 1
    draw_mc_scatter_panel(fig, M, lang, ...
        Ip_c, R_c, Ip_h, R_h, brine);
    add_mc_figure_legend(fig, lang, M, ...
        'scatter', brine);
    add_mc_condition_header(fig, lang, cn_note, en_note, []);
else
    Pe = draw_mc_lda_panel(fig, M, lang, Ip_c, R_c, Ip_h, R_h);
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
pdf_ci = calc_1d_kde(Ip_c, xi);
pdf_hi = calc_1d_kde(Ip_h, xi);
fill(ax_top, xi, pdf_ci, .24*M.col_ch4+.76*[1 1 1], ...
    'EdgeColor', M.col_ch4, 'LineWidth', 1.9, ...
    'HandleVisibility', 'off');
fill(ax_top, xi, pdf_hi, .24*M.col_h2+.76*[1 1 1], ...
    'EdgeColor', M.col_h2, 'LineWidth', 1.9, ...
    'HandleVisibility', 'off');
xlim(ax_top, xl);
ylim(ax_top, [0 1.06*max([pdf_ci(:); pdf_hi(:)])]);
axis(ax_top, 'off');

ax_side = axes('Parent', fig, 'Units', 'normalized', ...
    'Position', [main_x+main_w+side_gap main_y side_w main_h], ...
    'Color', 'white');
hold(ax_side, 'on');
yi = linspace(yl(1), yl(2), 220);
pdf_cr = calc_1d_kde(R_c, yi);
pdf_hr = calc_1d_kde(R_h, yi);
fill(ax_side, pdf_cr, yi, .24*M.col_ch4+.76*[1 1 1], ...
    'EdgeColor', M.col_ch4, 'LineWidth', 1.9, ...
    'HandleVisibility', 'off');
fill(ax_side, pdf_hr, yi, .24*M.col_h2+.76*[1 1 1], ...
    'EdgeColor', M.col_h2, 'LineWidth', 1.9, ...
    'HandleVisibility', 'off');
xlim(ax_side, [0 1.06*max([pdf_cr(:); pdf_hr(:)])]);
ylim(ax_side, yl);
axis(ax_side, 'off');
end

function Pe = draw_mc_lda_panel(fig, M, lang, Ip_c, R_c, Ip_h, R_h)
Xc = [Ip_c R_c]; Xh = [Ip_h R_h];
mu_c = mean(Xc, 1); mu_h = mean(Xh, 1);
Sw_cov = cov(Xc) + cov(Xh) + 1e-10*eye(2);
lda_direction = Sw_cov \ (mu_h - mu_c)';
z_c0 = Xc*lda_direction; z_h0 = Xh*lda_direction;
p_mean = .5*(mean(z_c0) + mean(z_h0));
p_std = max(sqrt(.5*(var(z_c0) + var(z_h0))), eps);
z_c = (z_c0-p_mean)/p_std;
z_h = (z_h0-p_mean)/p_std;
z_axis = linspace(min([z_c; z_h])-.8, max([z_c; z_h])+.8, 340);
pdf_c = calc_1d_kde(z_c, z_axis);
pdf_h = calc_1d_kde(z_h, z_axis);
pdf_ov = min(pdf_c, pdf_h);
Pe = .5*trapz(z_axis, pdf_ov);

ax = axes('Parent', fig, 'Units', 'normalized', ...
    'Position', [.145 .185 .735 .594], 'Color', 'white');
hold(ax, 'on');
fill(ax, z_axis, pdf_c, .20*M.col_ch4+.80*[1 1 1], ...
    'EdgeColor', M.col_ch4, 'LineWidth', 2.25, ...
    'HandleVisibility', 'off');
fill(ax, z_axis, pdf_h, .20*M.col_h2+.80*[1 1 1], ...
    'EdgeColor', M.col_h2, 'LineWidth', 2.25, ...
    'HandleVisibility', 'off');
fill(ax, z_axis, pdf_ov, [.93 .75 .77], ...
    'EdgeColor', [.80 .18 .25], 'LineWidth', 1.3, ...
    'HandleVisibility', 'off');
xline(ax, 0, ':', 'LineWidth', 1.5, ...
    'Color', [.32 .36 .41], 'HandleVisibility', 'off');
style_axes(ax, 18.3);
apply_axis_labels(ax, lang, ...
    '标准化 LDA 投影得分 η', 'Standardized LDA projection score η', ...
    '概率密度', 'Probability density');
ylim(ax, [0 1.07*max([pdf_c(:); pdf_h(:)])]);
grid(ax, 'on');
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
        'LDA 重叠密度', 'LDA overlap density', .275);
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
    [pe_text, pe_interp] = scientific_label( ...
        sprintf('P_e (KDE) = %.2f%%', Pe*100));
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

function export_bilingual_figure(fig, out_dir, base_name)
% Export BOTH vector PDF and vector SVG at the TRUE physical page size.
% Portrait outputs: 2.90 x 3.60 inches; landscape: 6.40 x 3.60.
% All standalone PDF/SVG pages have height = 2/5 of a 16 x 9 PPT slide.
if ~isgraphics(fig, 'figure')
    error('RPT:InvalidFigure', 'Invalid figure object: %s', base_name);
end
pdf_file = fullfile(out_dir, [base_name '.pdf']);
svg_file = fullfile(out_dir, [base_name '.svg']);
drawnow;
% Printing can follow colorbar repositioning, so lock the final physical
% typographic hierarchy just before PDF/SVG output.
lock_final_panel_typography(fig);
drawnow;
assert_no_panel_titles(fig);
validate_figure_typography(fig);
assert_consistent_main_axis_fonts(fig);
paper_in = getappdata(fig, 'PhysicalExportSizeInches');
assert(~isempty(paper_in) && numel(paper_in)==2 && ...
    abs(paper_in(2)-3.6)<1e-8, ...
    'RPT:UnexpectedPageSize', 'Every output must be 3.60 inches tall.');
set(fig, 'Renderer', 'painters', 'Color', 'white', ...
    'PaperUnits', 'inches', 'PaperSize', paper_in, ...
    'PaperPositionMode', 'manual', ...
    'PaperPosition', [0 0 paper_in(1) paper_in(2)], ...
    'InvertHardcopy', 'off');
% -r350 is only a fall-back resolution for any objects MATLAB rasterizes.
print(fig, pdf_file, '-dpdf', '-painters', '-r350');
fprintf('Saved PDF (%.2f x %.2f in): %s\n', paper_in(1), paper_in(2), pdf_file);
export_svg_figure(fig, svg_file);
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

function export_svg_figure(fig, svg_file)
% print -dsvg respects the explicitly assigned PaperPosition/PaperSize.
% exportgraphics is only a fallback and may crop to content bounds.
try
    print(fig, svg_file, '-dsvg', '-painters');
catch printErr
    try
        exportgraphics(fig, svg_file, 'ContentType', 'vector', ...
            'BackgroundColor', 'white');
        warning('RPT:SVGMayBeCropped', ...
            'SVG used exportgraphics fallback; verify its slide aspect ratio.');
    catch exportErr
        combined = MException('RPT:SVGExportFailed', ...
            'Failed to export SVG: %s', svg_file);
        combined = addCause(combined, printErr);
        combined = addCause(combined, exportErr);
        throwAsCaller(combined);
    end
end
fprintf('Saved vector SVG: %s\n', svg_file);
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
V.phi_dense = linspace(0.05, 0.30, 60);
V.Sg_dense  = linspace(0, 1, 60);
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
    .10, .50, A.T_rpt, A.P_rpt, ...
    'Methane', A.salinity, A.K_solid, A.rho_solid);
[V.Ip_h_ref, V.R_h_ref] = calc_rpt_phi_sg(A.Kd_10, A.Gd_10, ...
    .10, .50, A.T_rpt, A.P_rpt, ...
    'Hydrogen', A.salinity, A.K_solid, A.rho_solid);

V.Sg_list = [.20 .50 1.00];
V.P_sweep = linspace(1e6, 35e6, 120);
V.T_sweep = linspace(40, 120, 100);
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