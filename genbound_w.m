
%% Trajectory synthesis

dt = 0.01;

% generate upper and lower bounds
%  [yl, yu, xl, xu ] = genBounds(x, y)
track_ref = load("TestTrack.mat");
originPoint = [287; -176];
endPoint = TestTrack.cline(:,end);

R = vecnorm((TestTrack.bl - originPoint));

% [a, b, c] = linearInterpolate(T, index);
%     a*x + b*y = c

%%%%%%%%% z value needs to be clarified %%%%%%%%%%%%%%
%%%%%%%%% this function is based on z = [x_i y_i psi_i ... ] 
function [g, h, dg, dh] = nonlcon(z)
%vehicle dynamic constants
    m = 1400;
    Nw = 2;
    f = 0.01;
    Iz = 2667;
    a = 1.35;
    b = 1.45;
    By = 0.27;
    Cy = 1.2;
    Dy = 0.7;
    Ey = -1.6;
    Shy = 0;
    Svy = 0;
    g = 9.806;
    
 %vehicle dynamic first calculation
    Fzf = b/(a+b)*m*g;
    Fzr = a/(a+b)*m*g;
    Fmax = 0.7*m*g;
    
    nBin = 100;    
    dt = 0.05;
    g= zeros(2*(nBin+1),1);
    dg = zeros(2*(nBin+1),size(z,1));

    for i = 1:nBin+1
        r = norm([z(6*i+1); z(6*i+3)] - originPoint);
        index = find(r > R);
        [al,bl,cl] = linearInterpolate(TestTrack.bl,index);
        [ar,br,cr] = linearInterpolate(TestTrack.br,index);
        % g matrix
        g(2*i-1) = -cl+al*z(6*i+1)+bl*z(6*i+3);
        g(2*i) = cr-ar*z(6*i+1)-br*z(6*i+3);
        % dg matrix
        dg(2*i-1,6*i-5) = al;
        dg(2*i-1,6*i-4) = bl;
        dg(2*i,6*i-5) = -ar;
        dg(2*i,6*i-4) = -br;
    end
  
    h = zeros(6*(nBin+1),1);
    h(1:6) = z(1:6);
    for i = 1:nBin
        % h matrix
        u = z(6*i-4);
        v = z(6*i-2);
        psi = z(6*i-1);
        r = z(6*i);
        Fx = z(6*(1+nBin)+2*i-1);
        theta_f = z(6*(1+nBin)+2*i);
         
        alphaf = theta_f - atan((v+a*r)/u);
        alphar = -atan((v-b*r)/u);
        phiyf = (1-Ey)*(alphaf+Shy)+Ey/By*atan(By*(alphaf+Shy));
        phiyr = (1-Ey)*(alphar+Shy)+Ey/By*atan(By*(alphar+Shy));
        Fyf = Fzf*Dy*sin(Cy*atan(By*phiyf))+Svy;
        Fyr = Fzr*Dy*sin(Cy*atan(By*phiyr))+Svy;

        Ftotal = sqrt((Nw*Fx)^2+Fyr^2);
        if Ftotal>Fmax
            Fx = Fmax/Ftotal*Fx;
            Fyr = Fmax/Ftotal*Fyr;
        end
        
        dz = [u*cos(psi)-v*sin(psi);...
            1/m*(-f*m*g+Nw*Fx-Fyf*sin(theta_f))+v*r;...
            u*sin(psi)+v*cos(psi);...
            1/m*(Fyf*cos(theta_f)+Fyr)-u*r;...
            r;...
            1/Iz*(a*Fyf*cos(theta_f)-b*Fyr)];
        h(6*i+1:6*i+6) = z(6*i+1:6*i+6)-z(6*i-5:6*i)-dt*dz;
    end
    %generate dg = [...] g = [...], h =[...] dh = [...]

end

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