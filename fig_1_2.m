function fig_methods_setup_contact_model()
% FIG_METHODS_SETUP_CONTACT_MODEL
% MATLAB R2016a compatible
%
% Produces two separate figures:
%   Figure 1: 3D views (2x2)
%   Figure 2: matching 2D x-y projections (2x2)
%
% Internal units:
%   length = cm
%   time   = s
%
% Ecological calibration:
%   dt        = 0.02 s
%   meanSpeed = 25 cm/s
%   mean step = meanSpeed * dt = 0.5 cm
%
% For the Ahmed et al. (2021) 3D step-length law,
%   E[l] = 4*sigmaStep/sqrt(2*pi),
% so sigmaStep = E[l] * sqrt(2*pi)/4.

clear; clc; close all;
rng(2);

%% ========================= ECOLOGICAL PARAMETERS =========================
dt = 0.02;                 % s
meanSpeed = 25;            % cm/s
meanStep = meanSpeed * dt; % cm

% Convert mean step to sigma for the 3D chi step-length law
sigmaStep = meanStep * sqrt(2*pi) / 4;   % ~0.313 cm

maxSteps  = 300;           % total simulated duration = 6 s

%% ========================= GEOMETRY =========================
% Square sticky panel
W = 10;      % cm
H = 10;      % cm

% Initial release region
xStart     = -50;   % cm
releaseRad = 10;     % cm

%% ========================= MOVEMENT PARAMETERS =========================
kappa_p_crw  = 3.0;
kappa_b_crw  = 0.0;

kappa_p_bcrw = 3.0;
kappa_b_bcrw = 0.5;

%% ========================= CONTACT MODEL =========================
pStick     = 0.5;
kappa_post = 5.0;
epsShift   = 1e-3;   % cm

%% ========================= DISPLAY / DOMAIN =========================
nInitShow = 10;
nTrajShow = 10;

xEscapePos = 100;    % cm
xEscapeNeg = -100;   % cm
rEscape    = 100;    % cm

xPlotMin = -60;
xPlotMax =   10;
yzPlot   =  20;

%% ========================= PRECOMPUTE DATA =========================

% (a) initial points
pts = sample_release_points(nInitShow, xStart, releaseRad);

% (b) CRW trajectories (free movement only)
paths_crw = cell(nTrajShow,1);
for i = 1:nTrajShow
    path = simulate_one_trajectory_free( ...
        sigmaStep, maxSteps, xStart, releaseRad, ...
        kappa_p_crw, kappa_b_crw, ...
        xEscapePos, xEscapeNeg, rEscape, false, ...
        W, H);
    paths_crw{i} = path;
end

% (c) BCRW trajectories; bias is removed after first crossing of x = 0
paths_bcrw = cell(nTrajShow,1);
for i = 1:nTrajShow
    path = simulate_one_trajectory_free( ...
        sigmaStep, maxSteps, xStart, releaseRad, ...
        kappa_p_bcrw, kappa_b_bcrw, ...
        xEscapePos, xEscapeNeg, rEscape, true, ...
        W, H);
    paths_bcrw{i} = path;
end

% (d) one captured path and one edge-leakage path
foundCap  = false;
foundLeak = false;
nTry = 0;

while ~(foundCap && foundLeak) && nTry < 30000
    nTry = nTry + 1;

    [status, path] = simulate_one_trajectory_with_status_contact( ...
        sigmaStep, maxSteps, xStart, releaseRad, ...
        W, H, kappa_p_bcrw, kappa_b_bcrw, ...
        xEscapePos, xEscapeNeg, rEscape, true, ...
        pStick, kappa_post, epsShift);

    if status == 1 && ~foundCap
        path_cap = path;
        foundCap = true;
    elseif status == 2 && ~foundLeak
        if max(path(:,1)) > 1
            path_leak = path;
            foundLeak = true;
        end
    end
end

%% ========================= FIGURE 1: 3D =========================
figure('Color','w','Position',[60 60 1200 950]);

% (a)
subplot(2,2,1); hold on;
title('(a) Initial configuration', 'FontName', 'Arial', 'FontWeight', 'normal');
draw_panel_3d(W,H);
draw_release_disc_3d(xStart, releaseRad);
plot3(pts(:,1), pts(:,2), pts(:,3), 'ko', ...
    'MarkerSize', 5, 'MarkerFaceColor', [0.2 0.2 0.2]);
format_axes_3d(xPlotMin,xPlotMax,yzPlot);

% (b)
subplot(2,2,2); hold on;
title('(b) Correlated random walk', 'FontName', 'Arial', 'FontWeight', 'normal');
draw_panel_3d(W,H);
draw_release_disc_3d(xStart, releaseRad);
for i = 1:nTrajShow
    path = paths_crw{i};
    plot3(path(:,1), path(:,2), path(:,3), '-', 'LineWidth', 1.1);
end
format_axes_3d(xPlotMin,xPlotMax,yzPlot);

% (c)
subplot(2,2,3); hold on;
title('(c) Biased correlated random walk', 'FontName', 'Arial', 'FontWeight', 'normal');
draw_panel_3d(W,H);
draw_release_disc_3d(xStart, releaseRad);
for i = 1:nTrajShow
    path = paths_bcrw{i};
    plot3(path(:,1), path(:,2), path(:,3), '-', 'LineWidth', 1.1);
end
format_axes_3d(xPlotMin,xPlotMax,yzPlot);

% (d)
subplot(2,2,4); hold on;
title('(d) Panel-interaction outcomes', 'FontName', 'Arial', 'FontWeight', 'normal');
draw_panel_3d(W,H);
draw_release_disc_3d(xStart, releaseRad);
if foundCap
    plot3(path_cap(:,1), path_cap(:,2), path_cap(:,3), 'b-', 'LineWidth', 1.1);
end
if foundLeak
    plot3(path_leak(:,1), path_leak(:,2), path_leak(:,3), 'k-', 'LineWidth', 1.1);
end
format_axes_3d(xPlotMin,xPlotMax,yzPlot);

%% ========================= FIGURE 2: 2D x-y =========================
figure('Color','w','Position',[100 80 1200 950]);

% (a)
subplot(2,2,1); hold on;
title('(a) Initial configuration', 'FontName', 'Arial', 'FontWeight', 'normal');
draw_panel_xy(W);
draw_release_segment_xy(xStart, releaseRad);
plot(pts(:,1), pts(:,2), 'ko', ...
    'MarkerSize', 5, 'MarkerFaceColor', [0.2 0.2 0.2]);
format_axes_xy(xPlotMin,xPlotMax,yzPlot);

% (b)
subplot(2,2,2); hold on;
title('(b) Correlated random walk', 'FontName', 'Arial', 'FontWeight', 'normal');
draw_panel_xy(W);
draw_release_segment_xy(xStart, releaseRad);
for i = 1:nTrajShow
    path = paths_crw{i};
    plot(path(:,1), path(:,2), '-', 'LineWidth', 1.1);
end
format_axes_xy(xPlotMin,xPlotMax,yzPlot);

% (c)
subplot(2,2,3); hold on;
title('(c) Biased correlated random walk', 'FontName', 'Arial', 'FontWeight', 'normal');
draw_panel_xy(W);
draw_release_segment_xy(xStart, releaseRad);
for i = 1:nTrajShow
    path = paths_bcrw{i};
    plot(path(:,1), path(:,2), '-', 'LineWidth', 1.1);
end
format_axes_xy(xPlotMin,xPlotMax,yzPlot);

% (d)
subplot(2,2,4); hold on;
title('(d) Panel-interaction outcomes', 'FontName', 'Arial', 'FontWeight', 'normal');
draw_panel_xy(W);
draw_release_segment_xy(xStart, releaseRad);
if foundCap
    plot(path_cap(:,1), path_cap(:,2), 'b-', 'LineWidth', 1.1);
end
if foundLeak
    plot(path_leak(:,1), path_leak(:,2), 'k-', 'LineWidth', 1.1);
end
format_axes_xy(xPlotMin,xPlotMax,yzPlot);

end

%% =====================================================================
function format_axes_3d(xPlotMin,xPlotMax,yzPlot)
set(gca, 'FontName', 'Arial');
xlabel('x (cm)', 'FontName', 'Arial');
ylabel('y (cm)', 'FontName', 'Arial');
zlabel('z (cm)', 'FontName', 'Arial');
axis equal;
grid on;
box on;
view(3);
xlim([xPlotMin xPlotMax]);
ylim([-yzPlot yzPlot]);
zlim([-yzPlot yzPlot]);
end

%% =====================================================================
function format_axes_xy(xPlotMin,xPlotMax,yzPlot)
set(gca, 'FontName', 'Arial');
xlabel('x (cm)', 'FontName', 'Arial');
ylabel('y (cm)', 'FontName', 'Arial');
axis equal;
grid on;
box on;
xlim([xPlotMin xPlotMax]);
ylim([-yzPlot yzPlot]);
end

%% =====================================================================
function pts = sample_release_points(N, xStart, releaseRad)
pts = zeros(N,3);
for i = 1:N
    theta = 2*pi*rand;
    rr    = releaseRad*sqrt(rand);
    y0    = rr*cos(theta);
    z0    = rr*sin(theta);
    pts(i,:) = [xStart, y0, z0];
end
end

%% =====================================================================
function path = simulate_one_trajectory_free( ...
    sigmaStep, maxSteps, xStart, releaseRad, ...
    kappa_p, kappa_b, ...
    xEscapePos, xEscapeNeg, rEscape, useRandomPanelTarget, ...
    W, H)

theta = 2*pi*rand;
rr    = releaseRad*sqrt(rand);
y0    = rr*cos(theta);
z0    = rr*sin(theta);
r = [xStart, y0, z0];

v_prev = rand_unit_vector();

if useRandomPanelTarget && kappa_b > 0
    yTar = -W/2 + W*rand;
    zTar = -H/2 + H*rand;
    pTar = [0, yTar, zTar];
else
    pTar = [0, 0, 0]; %#ok<NASGU>
end

path = zeros(maxSteps + 1, 3);
idx = 1;
path(idx,:) = r;

attractionOn = useRandomPanelTarget && (kappa_b > 0);

for t = 1:maxSteps

    if attractionOn
        if useRandomPanelTarget
            g = pTar - r;
        else
            g = [1, 0, 0];
        end
        ng = norm(g);
        if ng < 1e-14
            g = [1,0,0];
        else
            g = g / ng;
        end
        kb_eff = kappa_b;
    else
        g = [0,0,0];
        kb_eff = 0;
    end

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

    % Remove panel-oriented bias after the first interaction with, or
    % passage through, the plane x = 0, in accordance with the Methods.
    if attractionOn
        [~, crossedPlane, ~] = segment_hits_panel(r, r_new, W, H);
        if crossedPlane
            attractionOn = false;
        end
    end

    idx = idx + 1;
    path(idx,:) = r_new;

    if r_new(1) > xEscapePos || r_new(1) < xEscapeNeg || norm(r_new) > rEscape
        path = path(1:idx,:);
        return;
    end

    r = r_new;
    v_prev = v_new;
end

path = path(1:idx,:);
end

%% =====================================================================
function [captured, path] = simulate_one_trajectory_contact( ...
    sigmaStep, maxSteps, xStart, releaseRad, ...
    W, H, kappa_p, kappa_b, ...
    xEscapePos, xEscapeNeg, rEscape, useRandomPanelTarget, ...
    pStick, kappa_post, epsShift)

theta = 2*pi*rand;
rr    = releaseRad*sqrt(rand);
y0    = rr*cos(theta);
z0    = rr*sin(theta);
r = [xStart, y0, z0];

v_prev = rand_unit_vector();

if useRandomPanelTarget && kappa_b > 0
    yTar = -W/2 + W*rand;
    zTar = -H/2 + H*rand;
    pTar = [0, yTar, zTar];
else
    pTar = [0, 0, 0]; %#ok<NASGU>
end

attractionOn = useRandomPanelTarget && (kappa_b > 0);

captured = false;
path = zeros(3*maxSteps + 20, 3);
idx = 1;
path(idx,:) = r;

for t = 1:maxSteps

    if attractionOn
        g = pTar - r;
        ng = norm(g);
        if ng < 1e-14
            g = [1,0,0];
        else
            g = g / ng;
        end
        kb_eff = kappa_b;
    else
        g = [0,0,0];
        kb_eff = 0;
    end

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
        idx = idx + 1;
        path(idx,:) = cp;

        if rand < pStick
            captured = true;
            path = path(1:idx,:);
            return;
        else
            % One outward step after failed contact, then CRW
            v_post = sample_post_contact_heading(kappa_post);
            ell_post = sample_step_length_eq10(sigmaStep);

            r = cp + epsShift * v_post;
            r_new2 = r + ell_post * v_post;

            idx = idx + 1;
            path(idx,:) = r_new2;

            if r_new2(1) > xEscapePos || r_new2(1) < xEscapeNeg || norm(r_new2) > rEscape
                path = path(1:idx,:);
                return;
            end

            r = r_new2;
            v_prev = v_post;
            attractionOn = false;
            continue;
        end

    elseif crossedPlane
        % After passing panel outside bounds: switch to CRW
        idx = idx + 1;
        path(idx,:) = cp;

        ell2 = sample_step_length_eq10(sigmaStep);
        r_new2 = cp + ell2 * v_new;

        idx = idx + 1;
        path(idx,:) = r_new2;

        if r_new2(1) > xEscapePos || r_new2(1) < xEscapeNeg || norm(r_new2) > rEscape
            path = path(1:idx,:);
            return;
        end

        r = r_new2;
        v_prev = v_new;
        attractionOn = false;
        continue;
    else
        idx = idx + 1;
        path(idx,:) = r_new;
    end

    if r_new(1) > xEscapePos || r_new(1) < xEscapeNeg || norm(r_new) > rEscape
        path = path(1:idx,:);
        return;
    end

    r = r_new;
    v_prev = v_new;
end

path = path(1:idx,:);
end

%% =====================================================================
function [status, path] = simulate_one_trajectory_with_status_contact( ...
    sigmaStep, maxSteps, xStart, releaseRad, ...
    W, H, kappa_p, kappa_b, ...
    xEscapePos, xEscapeNeg, rEscape, useRandomPanelTarget, ...
    pStick, kappa_post, epsShift)

% status:
%   1 = captured on panel
%   2 = edge leakage
%   3 = no contact before escape

theta = 2*pi*rand;
rr    = releaseRad*sqrt(rand);
y0    = rr*cos(theta);
z0    = rr*sin(theta);
r = [xStart, y0, z0];

v_prev = rand_unit_vector();

if useRandomPanelTarget && kappa_b > 0
    yTar = -W/2 + W*rand;
    zTar = -H/2 + H*rand;
    pTar = [0, yTar, zTar];
else
    pTar = [0, 0, 0]; %#ok<NASGU>
end

attractionOn = useRandomPanelTarget && (kappa_b > 0);

path = zeros(3*maxSteps + 20, 3);
idx = 1;
path(idx,:) = r;

leaked = false;

for t = 1:maxSteps

    if attractionOn
        g = pTar - r;
        ng = norm(g);
        if ng < 1e-14
            g = [1,0,0];
        else
            g = g / ng;
        end
        kb_eff = kappa_b;
    else
        g = [0,0,0];
        kb_eff = 0;
    end

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
        idx = idx + 1;
        path(idx,:) = cp;

        if rand < pStick
            status = 1;
            path = path(1:idx,:);
            return;
        else
            % One outward step after failed contact, then CRW
            v_post = sample_post_contact_heading(kappa_post);
            ell_post = sample_step_length_eq10(sigmaStep);

            r = cp + epsShift * v_post;
            r_new2 = r + ell_post * v_post;

            idx = idx + 1;
            path(idx,:) = r_new2;

            if r_new2(1) > xEscapePos || r_new2(1) < xEscapeNeg || norm(r_new2) > rEscape
                status = 3;
                path = path(1:idx,:);
                return;
            end

            r = r_new2;
            v_prev = v_post;
            attractionOn = false;
            continue;
        end

    elseif crossedPlane
        leaked = true;
        idx = idx + 1;
        path(idx,:) = cp;

        % Continue beyond the plane, then switch to CRW
        ell2 = sample_step_length_eq10(sigmaStep);
        r_new2 = cp + ell2 * v_new;

        idx = idx + 1;
        path(idx,:) = r_new2;

        if r_new2(1) > xEscapePos || r_new2(1) < xEscapeNeg || norm(r_new2) > rEscape
            status = 2;
            path = path(1:idx,:);
            return;
        end

        r = r_new2;
        v_prev = v_new;
        attractionOn = false;
        continue;
    else
        idx = idx + 1;
        path(idx,:) = r_new;
    end

    if r_new(1) > xEscapePos || r_new(1) < xEscapeNeg || norm(r_new) > rEscape
        if leaked
            status = 2;
        else
            status = 3;
        end
        path = path(1:idx,:);
        return;
    end

    r = r_new;
    v_prev = v_new;
end

if leaked
    status = 2;
else
    status = 3;
end
path = path(1:idx,:);
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
w = min(1,max(-1,w));

phi = 2*pi*rand;
s = sqrt(max(0,1-w^2));

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

v = cross(a,b);
s = norm(v);
c = dot(a,b);

if s < 1e-14
    if c > 0
        R = eye(3);
    else
        if abs(a(1)) < 0.9
            tmp = [1;0;0];
        else
            tmp = [0;1;0];
        end
        v = cross(a,tmp);
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

%% =====================================================================
function draw_panel_3d(W,H)
Y = [-W/2, W/2, W/2, -W/2, -W/2];
Z = [-H/2, -H/2, H/2, H/2, -H/2];
X = zeros(size(Y));

plot3(X,Y,Z,'r-','LineWidth',2);
fill3(X(1:4),Y(1:4),Z(1:4),[1.0 0.82 0.82], ...
    'FaceAlpha',0.35,'EdgeColor','none');
end

%% =====================================================================
function draw_release_disc_3d(xStart, releaseRad)
th = linspace(0,2*pi,200);
Y = releaseRad*cos(th);
Z = releaseRad*sin(th);
X = xStart*ones(size(th));

plot3(X,Y,Z,'b-','LineWidth',2);
fill3(X,Y,Z,[0.8 0.9 1.0], ...
    'FaceAlpha',0.25,'EdgeColor','none');
end

%% =====================================================================
function draw_panel_xy(W)
plot([0 0],[-W/2 W/2],'r-','LineWidth',3);
end

%% =====================================================================
function draw_release_segment_xy(xStart, releaseRad)
plot([xStart xStart],[-releaseRad releaseRad],'b-','LineWidth',3);
end