function dx = vehicle_model(x,u_in)

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

%data from x
u = x(2);
v = x(4);
psi = x(5);
r = x(6);

%data from u_in
delta_f = u_in(1);
Fx = u_in(2);

%vehicle dynamic first level calculation
Fzf = b/(a+b)*m*g;
Fzr = a/(a+b)*m*g;
Fmax = 0.7*m*g;

%second level calculation
alphaf = delta_f - atan((v+a*r)/u);
alphar = -atan((v-b*r)/u);
phiyf = (1-Ey)*(alphaf+Shy)+Ey/By*atan(By*(alphaf+Shy));
phiyr = (1-Ey)*(alphar+Shy)+Ey/By*atan(By*(alphar+Shy));
Fyf = Fzf*Dy*sin(Cy*atan(By*phiyf))+Svy;
Fyr = Fzr*Dy*sin(Cy*atan(By*phiyr))+Svy;

%Fx saturation
Ftotal = sqrt((Nw*Fx)^2+Fyr^2);
if Ftotal>Fmax
    Fx = Fmax/Ftotal*Fx;
    Fyr = Fmax/Ftotal*Fyr;
end

%dx calculation
dx = [u*cos(psi)-v*sin(psi);...
    1/m*(-f*m*g+Nw*Fx-Fyf*sin(delta_f))+v*r;...
    u*sin(psi)+v*cos(psi);...
    1/m*(Fyf*cos(delta_f)+Fyr)-u*r;...
    r;...
    1/Iz*(a*Fyf*cos(delta_f)-b*Fyr)];