function fig4_contact_frequency_2x2()
% FIG4_CONTACT_FREQUENCY_2X2
% Figure 4:
% Contact frequency and repeated panel interactions
%
% Top row:    pStick = 0
% Bottom row: pStick = 0.5
%
% Left panels:
%   Mean number of contacts per trajectory vs panel size L
%   for kappa_b = 0, 0.1, 0.2, 0.5, 1.0
%
% Right panels:
%   Unconditional distribution P(N = k) of contact counts
%   for fixed L = 10 and kappa_b = 0, 0.1, 0.2, 0.5, 1.0

clear; clc; close all;
rng(1);

%% ========================= PARAMETERS =========================
N = 1000;                  % use 100 for quick testing, 1000 for final

% Movement parameters (same baseline as Figure 3)
meanSpeed = 25;            % cm/s
dt        = 0.02;          % s
meanStep  = meanSpeed * dt;
sigmaStep = meanStep * sqrt(2*pi) / 4;

maxSteps  = 2500;
kappa_p   = 3.0;

% Bias values
kappa_b_vals_left  = [0.0, 0.1, 0.2, 0.5, 1.0];
kappa_b_vals_right = [0.0, 0.1, 0.2, 0.5, 1.0];

% Geometry
xStart     = -30;
releaseRad = 20;

% Panel sizes for left panels
L_vals = 2:2:30;

% Fixed panel size for right panels
L_fixed = 10;

% Post-contact parameters
kappa_post = 5.0;
epsShift   = 1e-3;

% Escape domain
xEscapePos = 100;
xEscapeNeg = -100;
rEscape    = 100;

% Stick probabilities for rows
pStick_rows = [0.0, 0.5];

%% ========================= STORAGE =========================
nRows    = length(pStick_rows);
nB_left  = length(kappa_b_vals_left);
nB_right = length(kappa_b_vals_right);
nL       = length(L_vals);

meanContacts = zeros(nRows, nB_left, nL);
distCounts   = cell(nRows, nB_right);

%% ========================= MAIN COMPUTATION =========================
for ir = 1:nRows
    pStick = pStick_rows(ir);

    % ----- Left panels: mean contacts vs L -----
    for ib = 1:nB_left
        kappa_b = kappa_b_vals_left(ib);

        for iL = 1:nL
            L = L_vals(iL);
            W = L;
            H = L;

            counts = zeros(N,1);

            for n = 1:N
                useRandomPanelTarget = (kappa_b > 0);

                counts(n) = simulate_one_trajectory_contact_count( ...
                    sigmaStep, maxSteps, xStart, releaseRad, ...
                    W, H, kappa_p, kappa_b, ...
                    xEscapePos, xEscapeNeg, rEscape, useRandomPanelTarget, ...
                    pStick, kappa_post, epsShift);
            end

            meanContacts(ir, ib, iL) = mean(counts);

            fprintf('Row %d: done left panel, pStick = %.1f, kappa_b = %.1f, L = %d\n', ...
                ir, pStick, kappa_b, L);
        end
    end

    % ----- Right panels: unconditional distribution at fixed L -----
    for ib = 1:nB_right
        kappa_b = kappa_b_vals_right(ib);
        W = L_fixed;
        H = L_fixed;

        counts = zeros(N,1);

        for n = 1:N
            useRandomPanelTarget = (kappa_b > 0);

            counts(n) = simulate_one_trajectory_contact_count( ...
                sigmaStep, maxSteps, xStart, releaseRad, ...
                W, H, kappa_p, kappa_b, ...
                xEscapePos, xEscapeNeg, rEscape, useRandomPanelTarget, ...
                pStick, kappa_post, epsShift);
        end

        distCounts{ir, ib} = counts;

        fprintf('Row %d: done right panel, pStick = %.1f, kappa_b = %.1f, L = %d\n', ...
            ir, pStick, kappa_b, L_fixed);
    end
end

%% ========================= PLOT =========================
figure('Color','w','Position',[80 80 1400 900]);

maxCountShown = 8;     % clip larger counts into final bin
kvals = 0:maxCountShown;

% Same visual style as Figure 3
lineColors = [ ...
    0.00 0.00 0.00; ...
    0.56 0.72 0.88; ...
    0.26 0.49 0.72; ...
    0.73 0.58 0.80; ...
    0.53 0.25 0.67];

lineStyles = {'-', '--', '-.', ':', '-'};

% -------------------- (a) pStick = 0, mean contacts vs L --------------------
subplot(2,2,1); hold on;

for ib = 1:nB_left
    plot(L_vals, squeeze(meanContacts(1,ib,:)), ...
        'Color', lineColors(ib,:), ...
        'LineStyle', lineStyles{ib}, ...
        'LineWidth', 2.0);
end

ylabel('mean number of contacts');
title('(a) p_{stick} = 0', 'FontWeight', 'normal');

xlim([2 30]);

set(gca, ...
    'FontName', 'Arial', ...
    'FontSize', 14, ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'LineWidth', 1.5, ...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off');

% -------------------- (b) pStick = 0, unconditional PMF --------------------
subplot(2,2,2); hold on;

for ib = 1:nB_right
    counts = distCounts{1,ib};
    pmf = compute_unconditional_pmf(counts, maxCountShown);

    plot(kvals, pmf, ...
        'Color', lineColors(ib,:), ...
        'LineStyle', lineStyles{ib}, ...
        'Marker', 'o', ...
        'LineWidth', 2.0, ...
        'MarkerSize', 5);
end

ylabel('probability');
title('(b) p_{stick} = 0', 'FontWeight', 'normal');

xlim([0 maxCountShown]);

set(gca, ...
    'FontName', 'Arial', ...
    'FontSize', 14, ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'LineWidth', 1.5, ...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off');

% -------------------- (c) pStick = 0.5, mean contacts vs L --------------------
subplot(2,2,3); hold on;

for ib = 1:nB_left
    plot(L_vals, squeeze(meanContacts(2,ib,:)), ...
        'Color', lineColors(ib,:), ...
        'LineStyle', lineStyles{ib}, ...
        'LineWidth', 2.0);
end

xlabel('panel side length, L');
ylabel('mean number of contacts');
title('(c) p_{stick} = 0.5', 'FontWeight', 'normal');

xlim([2 30]);

set(gca, ...
    'FontName', 'Arial', ...
    'FontSize', 14, ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'LineWidth', 1.5, ...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off');

% -------------------- (d) pStick = 0.5, unconditional PMF --------------------
subplot(2,2,4); hold on;

for ib = 1:nB_right
    counts = distCounts{2,ib};
    pmf = compute_unconditional_pmf(counts, maxCountShown);

    plot(kvals, pmf, ...
        'Color', lineColors(ib,:), ...
        'LineStyle', lineStyles{ib}, ...
        'Marker', 'o', ...
        'LineWidth', 2.0, ...
        'MarkerSize', 5);
end

xlabel('number of contacts per trajectory');
ylabel('probability');
title('(d) p_{stick} = 0.5', 'FontWeight', 'normal');

hleg = legend(make_bias_legend(kappa_b_vals_right), ...
    'Location', 'northeast');
set(hleg, 'FontSize', 14);
legend boxoff;

xlim([0 maxCountShown]);

set(gca, ...
    'FontName', 'Arial', ...
    'FontSize', 14, ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'LineWidth', 1.5, ...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off');

end

%% =====================================================================
function pmf = compute_unconditional_pmf(counts, maxCountShown)

kvals = 0:maxCountShown;
countsClip = counts;
countsClip(countsClip > maxCountShown) = maxCountShown;

pmf = zeros(size(kvals));

for k = 1:length(kvals)
    pmf(k) = sum(countsClip == kvals(k)) / length(countsClip);
end

end

%% =====================================================================
function labels = make_bias_legend(kappa_b_vals)

labels = cell(length(kappa_b_vals),1);

for i = 1:length(kappa_b_vals)
    labels{i} = sprintf('\\kappa_b = %.1f', kappa_b_vals(i));
end

end

%% =====================================================================
function nContacts = simulate_one_trajectory_contact_count( ...
    sigmaStep, maxSteps, xStart, releaseRad, ...
    W, H, kappa_p, kappa_b, ...
    xEscapePos, xEscapeNeg, rEscape, useRandomPanelTarget, ...
    pStick, kappa_post, epsShift)

% Initial position: uniform in release disc
theta = 2*pi*rand;
rr    = releaseRad * sqrt(rand);
y0    = rr * cos(theta);
z0    = rr * sin(theta);
r = [xStart, y0, z0];

% Initial heading: uniform on sphere
v_prev = rand_unit_vector();

% Fixed target point on panel for biased case
if useRandomPanelTarget && kappa_b > 0
    yTar = -W/2 + W*rand;
    zTar = -H/2 + H*rand;
    pTar = [0, yTar, zTar];
else
    pTar = [0, 0, 0];
end

% Bias acts only until first interaction with panel plane
attractionOn = useRandomPanelTarget && (kappa_b > 0);

nContacts = 0;

for t = 1:maxSteps

    % Attraction direction
    if attractionOn
        g = pTar - r;
        ng = norm(g);

        if ng < 1e-14
            g = [1, 0, 0];
        else
            g = g / ng;
        end

        kb_eff = kappa_b;

    else
        g = [0, 0, 0];
        kb_eff = 0;
    end

    % vMF heading update
    m = kappa_p * v_prev + kb_eff * g;
    kmag = norm(m);

    if kmag < 1e-14
        mu = rand_unit_vector();
        kappa = 0;
    else
        mu = m / kmag;
        kappa = kmag;
    end

    v_new = sample_vMF_3d(mu, kappa);
    ell   = sample_step_length_eq10(sigmaStep);
    r_new = r + ell * v_new;

    [isHit, crossedPlane, cp] = segment_hits_panel(r, r_new, W, H);

    if isHit

        % Contact with panel
        nContacts = nContacts + 1;

        if rand < pStick
            return;

        else
            % Failed contact: restart just behind panel, moving away
            v_post   = sample_post_contact_heading(kappa_post);
            ell_post = sample_step_length_eq10(sigmaStep);

            r = cp + epsShift * v_post;
            r_new2 = r + ell_post * v_post;

            if r_new2(1) > xEscapePos || ...
               r_new2(1) < xEscapeNeg || ...
               norm(r_new2) > rEscape
                return;
            end

            r = r_new2;
            v_prev = v_post;
            attractionOn = false;
            continue;
        end

    elseif crossedPlane

        % Crossed plane outside panel bounds: continue as CRW
        ell2 = sample_step_length_eq10(sigmaStep);
        r_new2 = cp + ell2 * v_new;

        if r_new2(1) > xEscapePos || ...
           r_new2(1) < xEscapeNeg || ...
           norm(r_new2) > rEscape
            return;
        end

        r = r_new2;
        v_prev = v_new;
        attractionOn = false;
        continue;

    else

        % Standard free-flight step
        if r_new(1) > xEscapePos || ...
           r_new(1) < xEscapeNeg || ...
           norm(r_new) > rEscape
            return;
        end

        r = r_new;
        v_prev = v_new;
    end
end

end

%% =====================================================================
function ell = sample_step_length_eq10(sigmaStep)

ell = sigmaStep * sqrt(randn^2 + randn^2 + randn^2);

end

%% =====================================================================
function v = sample_post_contact_heading(kappa_post)

mu_post = [-1, 0, 0];

while true
    v = sample_vMF_3d(mu_post, kappa_post);

    if v(1) < 0
        return;
    end
end

end

%% =====================================================================
function [tf, crossedPlane, point] = segment_hits_panel(r1, r2, W, H)

x1 = r1(1);
x2 = r2(1);

crossedPlane = false;
point = [NaN NaN NaN];

if (x1 < 0 && x2 < 0) || (x1 > 0 && x2 > 0)
    tf = false;
    return;
end

dx = x2 - x1;

if abs(dx) < 1e-14
    tf = false;
    return;
end

tau = -x1 / dx;

if tau < 0 || tau > 1
    tf = false;
    return;
end

crossedPlane = true;

yint = r1(2) + tau * (r2(2) - r1(2));
zint = r1(3) + tau * (r2(3) - r1(3));

point = [0, yint, zint];

tf = (abs(yint) <= W/2) && (abs(zint) <= H/2);

end

%% =====================================================================
function x = sample_vMF_3d(mu, kappa)

mu = mu(:)' / norm(mu);

if kappa < 1e-10
    x = rand_unit_vector();
    return;
end

u = rand;
a = exp(-kappa);
b = exp(kappa);

w = (1/kappa) * log(a + u * (b - a));
w = min(1, max(-1, w));

phi = 2*pi*rand;
s = sqrt(max(0, 1 - w^2));

xp = [s*cos(phi), s*sin(phi), w];
north = [0 0 1];

if norm(mu - north) < 1e-12

    x = xp;

elseif norm(mu + north) < 1e-12

    x = [xp(1), xp(2), -xp(3)];

else

    R = rotation_matrix_from_a_to_b(north, mu);
    x = (R * xp(:))';
    x = x / norm(x);

end

end

%% =====================================================================
function R = rotation_matrix_from_a_to_b(a, b)

a = a(:) / norm(a);
b = b(:) / norm(b);

v = cross(a, b);
s = norm(v);
c = dot(a, b);

if s < 1e-14

    if c > 0

        R = eye(3);

    else

        if abs(a(1)) < 0.9
            tmp = [1;0;0];
        else
            tmp = [0;1;0];
        end

        v = cross(a, tmp);
        v = v / norm(v);

        K = [ ...
             0    -v(3)  v(2); ...
             v(3)  0    -v(1); ...
            -v(2)  v(1)  0];

        R = eye(3) + 2*K*K;

    end

    return;
end

K = [ ...
     0    -v(3)  v(2); ...
     v(3)  0    -v(1); ...
    -v(2)  v(1)  0];

R = eye(3) + K + K*K * ((1 - c) / (s^2));

end

%% =====================================================================
function v = rand_unit_vector()

x = randn(1,3);
v = x / norm(x);

end