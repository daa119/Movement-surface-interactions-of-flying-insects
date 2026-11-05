function fig3_encounter_capture_vs_L_over_lp()

% FIG3_ENCOUNTER_CAPTURE_VS_L_OVER_LP
% Updated to use movement and contact parameters consistent with Figures 1 and 2
% Produces one figure with two panels:
%   (a) Encounter probability vs L / ell_p
%   (b) Capture probability vs L / ell_p

clear; clc; close all;
rng(1);

%% ========================= PARAMETERS =========================
N = 1000;                    % number of trajectories per panel size

% Ecological movement parameters (from Figures 1 & 2)
meanSpeed = 25;            % cm/s
dt        = 0.02;          % s
meanStep  = meanSpeed * dt; % cm
sigmaStep = meanStep * sqrt(2*pi)/4;  % step-length scale parameter

maxSteps  = 2500;          % total simulated steps

kappa_p   = 3.0;           % persistence parameter
kappa_b_vals = [0.0, 0.1, 0.2, 0.5, 1.0]; % bias strength values

% Geometry
xStart     = -30;          % initial release plane
releaseRad = 20;           % radius of circular release region

% Panel sizes
L_vals = 2:2:30;            % characteristic panel size (square panels)
% W = H = L in loop

% Contact model parameters
pStick     = 0.5;          % probability of capture
kappa_post = 5.0;          % post-contact concentration
epsShift   = 1e-3;         % small displacement away from panel

% Escape domain
xEscapePos = 100;
xEscapeNeg = -100;
rEscape    = 100;

%% ========================= PERSISTENCE LENGTH =========================
% Exact mean step length for 3D chi-step distribution
meanStep_exact = 4 * sigmaStep / sqrt(2*pi);

% Directional correlation for 3D vMF
rho = coth(kappa_p) - 1 / kappa_p;
ell_p = -meanStep_exact / log(rho);

fprintf('Persistence length ell_p = %.6f\n', ell_p);

%% ========================= MAIN LOOP =========================
encounter_prob = zeros(length(kappa_b_vals), length(L_vals));
capture_prob   = zeros(length(kappa_b_vals), length(L_vals));

for ib = 1:length(kappa_b_vals)

    kappa_b = kappa_b_vals(ib);

    for iL = 1:length(L_vals)

        L = L_vals(iL);
        W = L;
        H = L;

        nEncounter = 0;
        nCapture   = 0;

        for n = 1:N

            useRandomPanelTarget = (kappa_b > 0);

            [encFlag, capFlag] = simulate_one_trajectory_encounter_capture( ...
                sigmaStep, maxSteps, xStart, releaseRad, ...
                W, H, kappa_p, kappa_b, ...
                xEscapePos, xEscapeNeg, rEscape, useRandomPanelTarget, ...
                pStick, kappa_post, epsShift);

            if encFlag
                nEncounter = nEncounter + 1;
            end

            if capFlag
                nCapture = nCapture + 1;
            end

        end

        encounter_prob(ib, iL) = nEncounter / N;
        capture_prob(ib, iL)   = nCapture   / N;

        fprintf('Done: kappa_b = %.2f, L = %.1f\n', kappa_b, L);
    end
end

%% ========================= PLOT =========================
figure('Color','w','Position',[100 100 1100 500]);

% Colours for kappa_b values
line_colors = [
    0.00 0.00 0.00
    0.56 0.72 0.88
    0.26 0.49 0.72
    0.73 0.58 0.80
    0.53 0.25 0.67
];

% Line styles for kappa_b values
lineStyles = {'-', '--', '-.', ':', '-'};

hLeg = gobjects(length(kappa_b_vals),1);

% ---------- (a) Encounter probability ----------
subplot(1,2,1);
hold on;
axis square;

for ib = 1:length(kappa_b_vals)

    hLeg(ib) = plot(L_vals / ell_p, encounter_prob(ib,:), ...
        'LineStyle', lineStyles{ib}, ...
        'LineWidth', 2.0, ...
        'Color', line_colors(ib,:));

end

xlabel('panel size / persistence length, \nu', ...
    'FontSize', 12, ...
    'FontName', 'Arial', ...
    'Interpreter', 'tex');

ylabel('encounter probability', ...
    'FontSize', 12, ...
    'FontName', 'Arial', ...
    'Interpreter', 'tex');

%title('(a) Encounter probability', ...
%    'FontWeight', 'normal', ...
%    'FontSize', 12, ...
%    'FontName', 'Arial');

ylim([0 1]);
xlim([0, 25]);

grid off;

set(gca, ...
    'FontName', 'Arial', ...
    'FontSize', 14, ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'LineWidth', 1.5, ...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off');

hold off;


% ---------- (b) Capture probability ----------
subplot(1,2,2);
hold on;
axis square;

for ib = 1:length(kappa_b_vals)

    plot(L_vals / ell_p, capture_prob(ib,:), ...
        'LineStyle', lineStyles{ib}, ...
        'LineWidth', 2.0, ...
        'Color', line_colors(ib,:), ...
        'HandleVisibility', 'off');

end

xlabel('panel size / persistence length, \nu', ...
    'FontSize', 12, ...
    'FontName', 'Arial', ...
    'Interpreter', 'tex');

ylabel('capture probability', ...
    'FontSize', 12, ...
    'FontName', 'Arial', ...
    'Interpreter', 'tex');

%title('(b) Capture probability', ...
%    'FontWeight', 'normal', ...
%    'FontSize', 12, ...
%    'FontName', 'Arial');

ylim([0 1]);
xlim([0, 25]);

grid off;

set(gca, ...
    'FontName', 'Arial', ...
    'FontSize', 14, ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'LineWidth', 1.5, ...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off');
hold off;


%% Legend
leg_entries = cell(length(kappa_b_vals),1);

for ib = 1:length(kappa_b_vals)
    leg_entries{ib} = sprintf('\\kappa_b = %.2f', kappa_b_vals(ib));
end

lgd = legend(hLeg, leg_entries, ...
    'Interpreter', 'tex', ...
    'FontName', 'Arial', ...
    'FontSize', 14, ...
    'Box', 'off', ...
    'Location', 'best');

fprintf('\nSimulation complete.\n');

end


%% =====================================================================
function [encountered, captured] = simulate_one_trajectory_encounter_capture( ...
    sigmaStep, maxSteps, xStart, releaseRad, ...
    W, H, kappa_p, kappa_b, ...
    xEscapePos, xEscapeNeg, rEscape, useRandomPanelTarget, ...
    pStick, kappa_post, epsShift)

% Initial position
theta = 2*pi*rand;
rr    = releaseRad*sqrt(rand);
y0    = rr*cos(theta);
z0    = rr*sin(theta);
r = [xStart, y0, z0];

% Initial heading
v_prev = rand_unit_vector();

% Optional target point on the panel (fixed for entire trajectory)
if useRandomPanelTarget && kappa_b > 0
    yTar = -W/2 + W*rand;
    zTar = -H/2 + H*rand;
    pTar = [0, yTar, zTar];
else
    pTar = [0, 0, 0]; %#ok<NASGU>
end

% Attraction active only during the approach phase
attractionOn = useRandomPanelTarget && (kappa_b > 0);

encountered = false;
captured    = false;

for t = 1:maxSteps

    % Directional target during free flight
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

    % Heading update
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

        encountered = true;

        if rand < pStick

            captured = true;
            return;

        else

            % Failed contact: one outward step, then CRW
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

        % Passed the panel outside bounds: continue, but switch to CRW
        encountered = false;

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

yint = r1(2) + tau*(r2(2)-r1(2));
zint = r1(3) + tau*(r2(3)-r1(3));

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

w = (1/kappa) * log(a + u*(b-a));
w = min(1, max(-1, w));

phi = 2*pi*rand;
s = sqrt(max(0, 1-w^2));

xp = [s*cos(phi), s*sin(phi), w];
north = [0 0 1];

if norm(mu - north) < 1e-12

    x = xp;

elseif norm(mu + north) < 1e-12

    x = [xp(1), xp(2), -xp(3)];

else

    R = rotation_matrix_from_a_to_b(north, mu);
    x = (R*xp(:))';
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

        K = [   0   -v(3)  v(2);
              v(3)   0    -v(1);
             -v(2)  v(1)   0   ];

        R = eye(3) + 2*K*K;

    end

    return;
end

K = [   0   -v(3)  v(2);
      v(3)   0    -v(1);
     -v(2)  v(1)   0   ];

R = eye(3) + K + K*K*((1-c)/(s^2));

end


%% =====================================================================
function v = rand_unit_vector()

x = randn(1,3);
v = x / norm(x);

end