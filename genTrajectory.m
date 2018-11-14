%% Trajectory synthesis

dt = 0.01;

% generate upper and lower bounds
%  [yl, yu, xl, xu ] = genBounds(x, y)
track_ref = load("TestTrack.mat");
origin = [200; -200];

R = vecnorm((TestTrack.bl - origin));

[a, b, c] = linearInterpolate(T, index);
    a*x + b*y = c

function [g, h, dg, dh] = nonlcon(z)
    nBin = 100;
    i = 0:99;
    xIndex = 5*i + 1;
    for i = xIndex
        r = norm([z(i); z(i+2)] - origin);
        index = find(r > R);
        TestTrack.bl[index(end)] :TestTrack.bl[index(end) + 1]



%generate dg = [...] g = [...], h =[...] dh = [...]

%generate cost function J = [...]  dJ = [...]


%% MPC control
% linearize the system based on Task 1

% A(i), B(i) such that x(i+1) = A(i)*x(i) + B(i)*u(i)

% solve MPC controller (set Q R)
