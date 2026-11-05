function fig6_capture_vs_pstick()

clear; clc; close all;
rng(3);

%% ========================= PARAMETERS =========================
N = 1000;                 % trajectories per parameter combination

% Movement parameters (same baseline as Fig. 5)
meanSpeed = 25;           % cm/s
dt        = 0.02;         % s
meanStep  = meanSpeed * dt;
sigmaStep = meanStep * sqrt(2*pi) / 4;

maxSteps  = 2500;

% Geometry
xStart     = -30;
releaseRad = 20;
L          = 10;
W          = L;
H          = L;

% Post-contact parameters
kappa_post = 5.0;
epsShift   = 1e-3;

% Escape domain
xEscapePos = 100;
xEscapeNeg = -100;
rEscape    = 100;

% p_stick sweep
pStick_vals = 0:0.1:1.0;

% Panel (a): vary bias at representative persistence
kappa_p_fixed = 3.0;
kappa_b_vals_a = [0.0, 0.1, 0.2, 0.5, 1.0];

% Panel (b): vary persistence at representative bias
kappa_b_fixed = 0.2;
kappa_p_vals_b = [0.0, 0.5, 1.0, 3.0, 5.0];

captureProb_a = zeros(length(kappa_b_vals_a), length(pStick_vals));
captureProb_b = zeros(length(kappa_p_vals_b), length(pStick_vals));

%% ========================= MAIN SIMULATION =========================
fprintf('Running Figure 6 simulations...\n');

% ---------- Panel (a): varying kappa_b ----------
for ib = 1:length(kappa_b_vals_a)
    kappa_b = kappa_b_vals_a(ib);

    fprintf('\nPanel (a): kappa_b = %.2f\n', kappa_b);

    for ip = 1:length(pStick_vals)
        pStick = pStick_vals(ip);
        nCaptured = 0;

        for n = 1:N
            useRandomPanelTarget = (kappa_b > 0);

            wasCaptured = simulate_one_trajectory_capture_indicator( ...
                sigmaStep, maxSteps, xStart, releaseRad, ...
                W, H, kappa_p_fixed, kappa_b, ...
                xEscapePos, xEscapeNeg, rEscape, useRandomPanelTarget, ...
                pStick, kappa_post, epsShift);

            nCaptured = nCaptured + wasCaptured;
        end

        captureProb_a(ib, ip) = 100 * nCaptured / N;
    end
end

% ---------- Panel (b): varying kappa_p ----------
for ik = 1:length(kappa_p_vals_b)
    kappa_p = kappa_p_vals_b(ik);

    fprintf('\nPanel (b): kappa_p = %.2f\n', kappa_p);

    for ip = 1:length(pStick_vals)
        pStick = pStick_vals(ip);
        nCaptured = 0;

        for n = 1:N
            useRandomPanelTarget = (kappa_b_fixed > 0);

            wasCaptured = simulate_one_trajectory_capture_indicator( ...
                sigmaStep, maxSteps, xStart, releaseRad, ...
                W, H, kappa_p, kappa_b_fixed, ...
                xEscapePos, xEscapeNeg, rEscape, useRandomPanelTarget, ...
                pStick, kappa_post, epsShift);

            nCaptured = nCaptured + wasCaptured;
        end

        captureProb_b(ik, ip) = 100 * nCaptured / N;
    end
end

%% ========================= PLOT =========================
figure('Color','w','Position',[80 80 1200 450]);

% Colour scheme for kappa_b
lineColors_a = [ ...
    0.00 0.00 0.00; ...
    0.56 0.72 0.88; ...
    0.26 0.49 0.72; ...
    0.73 0.58 0.80; ...
    0.63 0.42 0.74];

% Different colour scheme for kappa_p
lineColors_b = [ ...
    0.00 0.00 0.00; ...
    0.76 0.88 0.70; ...
    0.50 0.72 0.46; ...
    0.28 0.58 0.30; ...
    0.28 0.58 0.30];

lineStyles = {'-', '--', '-.', ':', '-'};

% ---------- (a) ----------
subplot(1,2,1); hold on;

for ib = 1:length(kappa_b_vals_a)
    plot(pStick_vals, captureProb_a(ib,:), ...
        'Color', lineColors_a(ib,:), ...
        'LineStyle', lineStyles{ib}, ...
        'LineWidth', 2.0);
end

xlabel('sticking probability, p_{stick}');
ylabel('capture probability (%)');

title('(a) Varying bias strength \kappa_b at fixed \kappa_p = 3.0', ...
    'FontWeight', 'normal');

hleg = legend(arrayfun(@(x) sprintf('\\kappa_b = %.1f', x), ...
    kappa_b_vals_a, 'UniformOutput', false), ...
    'Location', 'southeast');
set(hleg, 'FontSize', 14);
legend boxoff;

set(gca, ...
    'FontName', 'Arial', ...
    'FontSize', 14, ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'LineWidth', 1.5, ...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off');

xlim([0 1]);
set(gca, 'XTick', 0:0.1:1);
ylim([0 100]);

% ---------- (b) ----------
subplot(1,2,2); hold on;

for ik = 1:length(kappa_p_vals_b)
    plot(pStick_vals, captureProb_b(ik,:), ...
        'Color', lineColors_b(ik,:), ...
        'LineStyle', lineStyles{ik}, ...
        'LineWidth', 2.0);
end

xlabel('sticking probability, p_{stick}');

title('(b) Varying persistence strength \kappa_p at fixed \kappa_b = 0.2', ...
    'FontWeight', 'normal');

hleg = legend(arrayfun(@(x) sprintf('\\kappa_p = %.1f', x), ...
    kappa_p_vals_b, 'UniformOutput', false), ...
    'Location', 'southeast');
set(hleg, 'FontSize', 14);
legend boxoff;

set(gca, ...
    'FontName', 'Arial', ...
    'FontSize', 14, ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'LineWidth', 1.5, ...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off');

xlim([0 1]);
set(gca, 'XTick', 0:0.1:1);
ylim([0 100]);

fprintf('\nSimulation complete.\n');

end

%% =====================================================================
function wasCaptured = simulate_one_trajectory_capture_indicator( ...
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

% Bias acts only until first interaction with the panel plane
attractionOn = useRandomPanelTarget && (kappa_b > 0);

wasCaptured = 0;

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
        if rand < pStick

            wasCaptured = 1;
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

        % Crossed x = 0 plane outside panel bounds: continue as CRW
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

% No crossing of x = 0 plane
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