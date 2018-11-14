%% Trajectory synthesis

dt = 0.01;

% generate upper and lower bounds
%  [yl, yu, xl, xu ] = genBounds(x, y)
track_ref = load("TestTrack.mat");
originPoint = [287; -176];
endPoint = TestTrack.cline(:,end);

R = vecnorm((TestTrack.bl - originPoint));

[a, b, c] = linearInterpolate(T, index);
    a*x + b*y = c

function [g, h, dg, dh] = nonlcon(z)
    nBin = 100;
    i = 0:99;
    xIndex = 5*i + 1;
    for i = xIndex
        r = norm([z(i); z(i+2)] - originPoint);
        index = find(r > R);
        TestTrack.bl[index(end)] :TestTrack.bl[index(end) + 1]



%generate dg = [...] g = [...], h =[...] dh = [...]

%generate cost function J = [...]  dJ = [...]
function [J, dJ] = costfun(z)
    nBin = 100;
    i = 0:(nBin-1);
    xIndex = 6*i + 1;
    uIndex = 6*i + 2;
    yIndex = 6*i + 3;
    vIndex = 6*i + 4;
    psiIndex = 6*i + 5;
    rIndex = 6*i + 6;
    end1 = 6*(nBin-1) + 6;    
    j = 0:(nBin-2);
    deltaIndex = end1 + 2*j + 1;
    FxIndex = end1 + 2*j + 2;
    
    x = zeros(size(i,2),1);
    u = zeros(size(i,2),1);
    y = zeros(size(i,2),1);
    v = zeros(size(i,2),1);
    psi = zeros(size(i,2),1);
    r = zeros(size(i,2),1);
    delta = zeros(size(j,2),1);
    Fx = zeros(size(j,2),1);
    
    x(i+1,1) = z(xIndex);
    u(i+1,1) = z(uIndex);
    y(i+1,1) = z(yIndex);
    v(i+1,1) = z(vIndex);
    psi(i+1,1) = z(psiIndex);
    delta(j+1,1) = z(deltaIndex);
    Fx(j+1,1) = z(FxIndex);
    
    % size of J must be 1 x 1
    J1 = zeros(size(i,2),1);
    J2 = zeros(size(j,2),1);
    
    J1(i+1,1) = (x(i+1,1)-endPoint(1,1)).^2 + u(i+1,1).^2 + (y(i+1,1)-endPoint(2,1)).^2 + v(i+1,1).^2 + psi(i+1,1).^2 + r(i+1,1).^2;
    J2(j+1,1) = delta(j+1,1).^2 + Fx(j+1,1).^2;
    
    J = sum(J1) + sum(J2);
    
    % size of dJ must be 1 x 603 (1 x no. of elements in 'z')
    dJ = 2.0*z';
    dJ(xIndex) = dJ(xIndex) - 2.*endPoint(1,1);
    dJ(yIndex) = dJ(yIndex) - 2.*endPoint(2,1);
    
end


%% MPC control
% linearize the system based on Task 1

% A(i), B(i) such that x(i+1) = A(i)*x(i) + B(i)*u(i)

% solve MPC controller (set Q R)
