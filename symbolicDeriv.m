clear
load('vehicle_parameter.mat');
syms X u Y v Psi r Fx deltaf Fyf(phiyf) Fyr(phiyr) phiyf(alphaf) phiyr(alphar) alphaf(deltaf, v, r, u) alphar(v, r, u);

alphaf = (deltaf - atan( (v+a*r)/u))*180.0/pi;
alphar = -atan((v-b*r)/u)*180.0/pi;
phiyf = (1 - Ey)*alphaf + Ey/By*atan(By*alphaf);
phiyr = (1 - Ey)*alphar + Ey/By*atan(By*alphar);
Fyf = Fzf*Dy*sin(Cy*atan(By*phiyf));
Fyr = Fzr*Dy*sin(Cy*atan(By*phiyr));

global dz;
dz  = [ ...
            u*cos(Psi) - v*sin(Psi);
            1.0/m*(-f*m*g + Nw*Fx - Fyf*sin(deltaf)) + v*r;
            u*sin(Psi) + v*cos(Psi);
            1.0/m*(Fyf*cos(deltaf) + Fyr) - u*r;
            r;
            1.0/Iz*(a*Fyf*cos(deltaf) - b*Fyr) ...
         ];

A = jacobian(dz, [X, u, Y, v, Psi, r]);
B = jacobian(dz, [deltaf, Fx]);

states = zeros(100, 6);
input = zeros(99, 2);
input(:,1) = -pi/10.0;

states(1,:) = [287,5, -176, 0, 2, 0];

dt = 0.01;
for i = 1:99
    states(i+1, :) = states(i, :) + dt*carModel(states(i, :), input(i,:))';
end

plot(states(:,1), states(:,3), '-x')

function deriv = carModel(state, input)
    global dz;
    deriv = subs(dz, symvar(dz), [input(2), state(5), input(1), state(6), state(2), state(4)]);
end

