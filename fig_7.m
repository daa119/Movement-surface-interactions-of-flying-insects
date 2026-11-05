function fig7_time_to_first_encounter()

clear; clc; close all;
rng(4);

%% ========================= PARAMETERS =========================
N = 1000;                 % trajectories per parameter combination

% Movement parameters (same baseline as Fig. 6)
meanSpeed = 25;           % cm/s
dt        = 0.02;         % s
meanStep  = meanSpeed * dt;
sigmaStep = meanStep * sqrt(2*pi) / 4;

maxSteps  = 2500;

% Geometry
releaseRad = 20;
L          = 10;
W          = L;
H          = L;

% Escape domain
xEscapePos = 100;
xEscapeNeg = -100;
rEscape    = 100;

% No retention here
pStick = 0.0; %#ok<NASGU>

% Initial offsets to test
xStart_vals = -10:-5:-100;    % cm
xOffset_vals = abs(xStart_vals);

% Panel (a): vary bias at representative persistence
kappa_p_fixed = 3.0;
kappa_b_vals_a = [0.0, 0.1, 0.2, 0.5, 1.0];

% Panel (b): vary persistence at representative bias
kappa_b_fixed = 0.2;
kappa_p_vals_b = [0.0, 0.5, 1.0, 3.0, 5.0];

meanTime_a = nan(length(kappa_b_vals_a), length(xStart_vals));
meanTime_b = nan(length(kappa_p_vals_b), length(xStart_vals));

%% ========================= MAIN SIMULATION =========================
fprintf('Running Figure 7 simulations...\n');

% ---------- Panel (a): varying kappa_b ----------
for ib = 1:length(kappa_b_vals_a)
    kappa_b = kappa_b_vals_a(ib);

    fprintf('\nPanel (a): kappa_b = %.2f\n', kappa_b);

    for ix = 1:length(xStart_vals)
        xStart = xStart_vals(ix);

        firstTimes = nan(N,1);

        for n = 1:N
            useRandomPanelTarget = (kappa_b > 0);

            tFirst = simulate_one_trajectory_first_encounter_time( ...
                sigmaStep, dt, maxSteps, xStart, releaseRad, ...
                W, H, kappa_p_fixed, kappa_b, ...
                xEscapePos, xEscapeNeg, rEscape, useRandomPanelTarget);

            firstTimes(n) = tFirst;
        end

        meanTime_a(ib, ix) = mean(firstTimes(~isnan(firstTimes)));
    end
end

% ---------- Panel (b): varying kappa_p ----------
for ik = 1:length(kappa_p_vals_b)
    kappa_p = kappa_p_vals_b(ik);

    fprintf('\nPanel (b): kappa_p = %.2f\n', kappa_p);

    for ix = 1:length(xStart_vals)
        xStart = xStart_vals(ix);

        firstTimes = nan(N,1);

        for n = 1:N
            useRandomPanelTarget = (kappa_b_fixed > 0);

            tFirst = simulate_one_trajectory_first_encounter_time( ...
                sigmaStep, dt, maxSteps, xStart, releaseRad, ...
                W, H, kappa_p, kappa_b_fixed, ...
                xEscapePos, xEscapeNeg, rEscape, useRandomPanelTarget);

            firstTimes(n) = tFirst;
        end

        meanTime_b(ik, ix) = mean(firstTimes(~isnan(firstTimes)));
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
    0.12 0.42 0.62];

lineStyles = {'-', '--', '-.', ':', '-'};

% ---------- (a) ----------
subplot(1,2,1); hold on;

for ib = 1:length(kappa_b_vals_a)
    plot(xOffset_vals, meanTime_a(ib,:), ...
        'Color', lineColors_a(ib,:), ...
        'LineStyle', lineStyles{ib}, ...
        'LineWidth', 2.0);
end

xlabel('initial distance from panel, {\it d} (cm)');
ylabel('mean time to first encounter, {\it t}_e (s)');

title('(a) Varying bias strength \kappa_b at fixed \kappa_p = 3.0', ...
    'FontWeight', 'normal');

hleg = legend(arrayfun(@(x) sprintf('\\kappa_b = %.1f', x), ...
    kappa_b_vals_a, 'UniformOutput', false), ...
    'Location', 'northwest');
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

xlim([min(xOffset_vals) max(xOffset_vals)]);

% ---------- (b) ----------
subplot(1,2,2); hold on;

for ik = 1:length(kappa_p_vals_b)
    plot(xOffset_vals, meanTime_b(ik,:), ...
        'Color', lineColors_b(ik,:), ...
        'LineStyle', lineStyles{ik}, ...
        'LineWidth', 2.0);
end

xlabel('initial distance from panel, {\it d} (cm)');

title('(b) Varying persistence strength \kappa_p at fixed \kappa_b = 0.2', ...
    'FontWeight', 'normal');

hleg = legend(arrayfun(@(x) sprintf('\\kappa_p = %.1f', x), ...
    kappa_p_vals_b, 'UniformOutput', false), ...
    'Location', 'northwest');
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

xlim([min(xOffset_vals) max(xOffset_vals)]);

fprintf('\nSimulation complete.\n');

end

%% =====================================================================
function tFirst = simulate_one_trajectory_first_encounter_time( ...
    sigmaStep, dt, maxSteps, xStart, releaseRad, ...
    W, H, kappa_p, kappa_b, ...
    xEscapePos, xEscapeNeg, rEscape, useRandomPanelTarget)

% Returns:
%   tFirst = time (in s) of first encounter with the finite panel
%          = NaN if no panel encounter occurs before leakage/escape

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

tFirst = NaN;

for t = 1:maxSteps

    % Attraction direction
    if kappa_b > 0
        if useRandomPanelTarget
            g = pTar - r;
        else
            g = [1, 0, 0];
        end

        ng = norm(g);

        if ng < 1e-14
            g = [1, 0, 0];
        else
            g = g / ng;
        end
    else
        g = [0, 0, 0];
    end

    % vMF heading update
    m = kappa_p * v_prev + kappa_b * g;
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

    [isHit, crossedPlane, ~] = segment_hits_panel(r, r_new, W, H);

    % Successful first encounter
    if isHit
        tFirst = t * dt;
        return;
    end

    % Edge leakage terminates trajectory without encounter
    if crossedPlane && ~isHit
        return;
    end

    % Escape without encounter
    if r_new(1) > xEscapePos || ...
       r_new(1) < xEscapeNeg || ...
       norm(r_new) > rEscape

        return;
    end

    r = r_new;
    v_prev = v_new;
end

end

%% =====================================================================
function ell = sample_step_length_eq10(sigmaStep)

ell = sigmaStep * sqrt(randn^2 + randn^2 + randn^2);

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