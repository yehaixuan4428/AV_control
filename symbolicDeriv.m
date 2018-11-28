clear
%load('vehicle_parameter.mat');

syms X u Y v Psi r Fx deltaf Fyf(phiyf) Fyr(phiyr) phiyf(alphaf) phiyr(alphar) alphaf(deltaf, v, r, u) alphar(v, r, u) m Nw f Iz a b By Cy Dy Ey g Fzf Fzr;

Fzf = b/(a+b)*m*g;
Fzr = a/(a+b)*m*g;
alphaf = deltaf - 1.0/tan( (v+a*r)/u);
alphar = -1.0/tan( (v-b*r)/u);
phiyf = (1 - Ey)*alphaf + Ey/By/tan(By*alphaf);
phiyr = (1 - Ey)*alphar + Ey/By/tan(By*alphar);
Fyf = Fzf*Dy*sin(Cy/tan(By*phiyf));
Fyr = Fzr*Dy*sin(Cy/tan(By*phiyr));

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
