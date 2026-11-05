function fig8_edge_leakage_heatmaps_kappap_kappab()
% FIG8_EDGE_LEAKAGE_HEATMAPS_KAPPAP_KAPPAB
% MATLAB R2016a compatible

clear; clc; close all;
rng(1);

%% ========================= PARAMETERS =========================
N = 1000;                 % trajectories per parameter combination

% Movement parameters
meanSpeed = 25;
dt        = 0.02;
meanStep  = meanSpeed * dt;
sigmaStep = meanStep * sqrt(2*pi) / 4;

maxSteps  = 2500;

% Geometry / release
xStart     = -30;
releaseRad = 20;

% Panel sizes for the three subpanels
L_vals = [10, 20, 30];
nL     = length(L_vals);

% Post-contact parameters
kappa_post = 5.0;
epsShift   = 1e-3;

% Escape domain
xEscapePos = 100;
xEscapeNeg = -100;
rEscape    = 100;

% Set retention to zero
pStick = 0.0; 

% Heat-map parameter grids
% Same parameter space as Figure 5
kappa_p_vals = 0:0.25:5.0;
kappa_b_vals = 0:0.025:0.5;

nKp = length(kappa_p_vals);
nKb = length(kappa_b_vals);

% Storage
leakProb = zeros(nKb, nKp, nL);

%% ========================= MAIN SIMULATION =========================
fprintf('Running edge-leakage heat maps...\n');

for iL = 1:nL

    L = L_vals(iL);
    W = L;
    H = L;

    fprintf('\n=== L = %d ===\n', L);

    for ikb = 1:nKb

        kappa_b = kappa_b_vals(ikb);

        for ikp = 1:nKp

            kappa_p = kappa_p_vals(ikp);

            nLeak = 0;

            for n = 1:N

                useRandomPanelTarget = (kappa_b > 0);

                wasLeak = simulate_one_trajectory_edge_leakage_indicator( ...
                    sigmaStep, maxSteps, xStart, releaseRad, ...
                    W, H, kappa_p, kappa_b, ...
                    xEscapePos, xEscapeNeg, rEscape, useRandomPanelTarget, ...
                    pStick, kappa_post, epsShift);

                nLeak = nLeak + wasLeak;

            end

            leakProb(ikb, ikp, iL) = nLeak / N;

        end

        fprintf('  Finished kappa_b = %.3f\n', kappa_b);

    end
end

%% ========================= PLOT =========================
figure('Color','w','Position',[60 80 1500 480]);

for iL = 1:nL

    subplot(1,3,iL);

    imagesc(kappa_p_vals, kappa_b_vals, ...
        100*leakProb(:,:,iL));

    set(gca,'YDir','normal');

    hold on;

    [C,h] = contour(kappa_p_vals, kappa_b_vals, ...
        100*leakProb(:,:,iL), ...
        [25 50 75], ...
        'k', ...
        'LineWidth',1.2);

    clabel(C,h, ...
        'FontSize',12, ...
        'Color','k');

    colormap(jet);
    caxis([0 100]);

    xlabel('\kappa_p (directional persistence)');

    if iL == 1
        ylabel('\kappa_b (bias strength)');
    end

    title(sprintf('L = %d cm',L_vals(iL)), ...
        'FontWeight','normal');

    set(gca, ...
        'FontName','Arial', ...
        'FontSize',14, ...
        'Box','off', ...
        'TickDir','out', ...
        'LineWidth',1.5, ...
        'XMinorTick','off', ...
        'YMinorTick','off');

    axis square;

    % Colour bar only on right-hand panel
    if iL == nL

        axPos = get(gca,'Position');

        c = colorbar;
        c.Label.String   = 'edge leakage probability (%)';
        c.Label.FontName = 'Arial';
        c.Label.FontSize = 14;

        % Restore third subplot to its original size
        set(gca,'Position',axPos);

        % Position colour bar to the right
        cPos = get(c,'Position');
        cPos(1) = axPos(1) + axPos(3) + 0.015;
        set(c,'Position',cPos);

    end

end

fprintf('\nSimulation complete.\n');

end

%% =====================================================================
function wasLeak = simulate_one_trajectory_edge_leakage_indicator( ...
    sigmaStep, maxSteps, xStart, releaseRad, ...
    W, H, kappa_p, kappa_b, ...
    xEscapePos, xEscapeNeg, rEscape, useRandomPanelTarget, ...
    pStick, kappa_post, epsShift) %#ok<INUSD>

% Returns:
%   wasLeak = 1 if the trajectory crosses the plane x = 0 outside the panel
%             bounds before any contact with the panel
%           = 0 otherwise

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

wasLeak = 0;

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

    ell = sample_step_length_eq10(sigmaStep);

    r_new = r + ell * v_new;

    [isHit, crossedPlane, ~] = ...
        segment_hits_panel(r, r_new, W, H);

    % Panel contact occurs before leakage
    if isHit
        return;
    end

    % Plane crossing outside the finite panel = edge leakage
    if crossedPlane
        wasLeak = 1;
        return;
    end

    % Escape without panel-plane crossing
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

if (x1 < 0 && x2 < 0) || ...
   (x1 > 0 && x2 > 0)

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

tf = (abs(yint) <= W/2) && ...
     (abs(zint) <= H/2);

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