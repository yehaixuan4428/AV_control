
%% parameters
%wheelbases
L=3;
%distance from rear wheel to center of mass
b=1.5;

%time discretization
dt=0.01;
%time span
T=0:dt:6;

%% load reference trajectory
load('part1_traj_05_timestep.mat')

U_ref=interp1(0:0.05:6,[U,U(:,end)]',T)';

Y_ref=interp1(0:0.05:6,Y,T)';

%the reference trajectory Yref is given in a 3x601 double

%% 2.1 Discrete-time A and B matrices
%these are the system linearized in discrete time about the reference
%trajectory i.e. x(i+1)=A_i*x_i+B_i*u_i

dxA=@(i) [
            0, 0, -U_ref(1,i)*sin(Y_ref(3,i))-b/L*U_ref(1,i)*tan(U_ref(2,i))*cos(Y_ref(3,i));
            0, 0,  U_ref(1,i)*cos(Y_ref(3,i))-b/L*U_ref(1,i)*tan(U_ref(2,i))*sin(Y_ref(3,i));
            0, 0, 0
       ];


dB=@(i) [ 
            cos(Y_ref(3,i))-b/L*tan(U_ref(2,i))*sin(Y_ref(3,i)), -b/L*U_ref(1,i)/(cos(U_ref(2,i)))^2*sin(Y_ref(3,i));
            sin(Y_ref(3,i))+b/L*tan(U_ref(2,i))*cos(Y_ref(3,i)), b/L*U_ref(1,i)/(cos(U_ref(2,i)))^2*cos(Y_ref(3,i));
            1.0/L*tan(U_ref(2,i)), U_ref(1,i)/(cos(U_ref(2,i)))^2/L
       ];

A = @(i) eye(3) + dt.*dxA(i);
B = @(i) dt.*dB(i);

%A(2)

%% 2.2 Number of decision variables for colocation method
npred=10;

Ndec= (npred + 1)*3 + npred*2;

%% 2.3 write and test function  to construct Aeq beq (equality constraints
%enforce x(i+1)=A_i*x_i+B_i*u_i for the prediction horizon) 
%check the values of Aeq and beq timestep 1 (t=0)
%with the given initial condition of Y(0)=[0.25;-0.25;-0.1];
[Aeq_test1 beq_test1] = eq_cons(1, A, B, npred, [0.25;-0.25;-0.1], 3, 2);
%Aeq_test1=
%beq_test1=

%% 2.4 write and test function to generate limits on inputs 
%check the value at timestep 1 (t=0)
[Lb_test1, Ub_test1] = bound_cons(1, U_ref, Y_ref, npred);

% 2.5 simulate controller working from initial condition [0.25;-0.25;-0.1]
%use ode45 to between inputs
Ysim = zeros(3, length(T)-npred+1);
Ytotal = zeros(3, length(T)-npred+1);
Ysim(:, 1) = [0.25;-0.25;-0.1];
Ytotal(:, 1) = Ysim(:,1) + Y_ref(:,1);
Q = [1 0 0; 0 1 0;0 0 0.5];
R = [0.1 0; 0 0.01];
nState = size(Q,1);
nInput = size(R,1);

j = 0:npred;
xIndex = 1 + nState*j;
yIndex = xIndex + 1;
psiIndex = xIndex + 2;
j = 0:npred-1;
uIndex = psiIndex(end) + 1 + nInput*j;
deltaIndex = uIndex + 1;

H = zeros(nState*(npred+1) + nInput*npred);
c = zeros(nState*(npred+1) + nInput*npred, 1);
for j = 1:npred
    H(xIndex(j):psiIndex(j), xIndex(j):psiIndex(j)) = Q;
    H(uIndex(j):deltaIndex(j), uIndex(j):deltaIndex(j)) = R;
end
H(xIndex(npred+1):psiIndex(npred+1), xIndex(npred+1):psiIndex(npred+1)) = Q;

for i = 1:length(T)-npred

    [Aeq, beq] = eq_cons(i, A, B, npred, Ysim(:, i), 3, 2);
    [Lb, Ub] = bound_cons(i, U_ref, Y_ref, npred);

    solDV = quadprog(H, c, [], [], Aeq, beq, Lb, Ub);
    %u = solDV( nState*(npred+1)+1:nState*(npred+1)+nInput);
    u = solDV(uIndex(1):deltaIndex(1));
    u = U_ref(:,i) + u;

    %if solDV(xIndex(1):psiIndex(1)) ~= Ysim(:,i)
        %error("first constraints wrong %d", i);
    %end

    %for k = 1:length(uIndex)
        %predict = A(i + k - 1)*solDV(xIndex(k):psiIndex(k)) + B(i + k - 1)*solDV(uIndex(k):deltaIndex(k));
        %if solDV(xIndex(k+1):psiIndex(k+1)) ~= predict
            %Aeq(xIndex(k+1):psiIndex(k+1), xIndex(k):deltaIndex(k+1))
            %beq(xIndex(k+1):psiIndex(k+1), :)
            %error("constraints wrong %d %d", i, k);
        %end
    %end

    [T, Y] = ode45(@(t,x)  odefun(x, u), [0,0.01], Ytotal(:,i));
    Ytotal(:,i+1) = Y(end,:)';
    Ysim(:,i+1) = Ytotal(:,i+1) - Y_ref(:,i+1);
end

%caclulate max distance error between the actual and nominal trajectories,
%when the x value of the actual trajectory is between or equal to 3 and 4 m
%index = find(Y_ref(1,1:length(T)-20) + Ysim(1,:) >= 3 & Y_ref(1,1:length(T)-20) + Ysim(1,:) <= 4);
index = find(Ytotal(1,:) >= 3 & Ytotal(1,:) <= 4);
errors = sqrt( Ysim(1,index).^2 + Ysim(2,index).^2);
[max_dist_error, maxID] = max(errors)

plot(Y_ref(1,:), Y_ref(2,:))
hold on
plot(Ytotal(1,:), Ytotal(2,:))
%hold on
%theta = 0:0.01:2*pi;
%hold on
%plot((0.7*cos(theta)+3.5),(0.7*sin(theta)-0.5))
legend('ref', 'total')


%% I highly recommend you write functions for constraint generation and dynamics down here and 
%call them when needed above, for exa

function [Aeq,beq]=eq_cons(idx,A,B,npred, x0, nState, nInput)
    i = 0:npred;
    xIndex = 1 + nState*i;
    yIndex = xIndex + 1;
    psiIndex = xIndex + 2;
    i = 0:npred-1;
    uIndex = psiIndex(end) + 1 + nInput*i;
    deltaIndex = uIndex + 1;

    Aeq = zeros(psiIndex(end), deltaIndex(end));
    beq = zeros(psiIndex(end), 1);
    Aeq(1:nState, 1:nState) = eye(nState);
    beq(1:nState) = x0;

    for i = 1:npred
        Aeq(xIndex(i+1):psiIndex(i+1), xIndex(i):psiIndex(i+1)) = [A(idx+i-1), -eye(3)];
        Aeq(xIndex(i+1):psiIndex(i+1), uIndex(i):deltaIndex(i)) = B(idx+i-1);
    end

%build matrix for A_i*x_i+B_i*u_i-x_{i+1}=0
%in the form Aeq*z=beq
%initial_idx specifies the time index of initial condition from the reference trajectory 
%A and B are function handles above
end

function [Lb,Ub]=bound_cons(idx,U_ref,Y_ref, npred)
%initial_idx is the index along uref the initial condition is at
    Ub = [repmat([inf;inf;inf], npred+1, 1); repmat([1;0.5], npred, 1)];
    Lb = [repmat([-inf;-inf;-inf], npred+1, 1); repmat([0;-0.5], npred, 1)];

    U_slice = reshape(U_ref(:, idx:idx+npred-1), [], 1);
    Y_slice = reshape(Y_ref(:, idx:idx+npred), [], 1);

    Ub = Ub - [Y_slice; U_slice];
    Lb = Lb - [Y_slice; U_slice];
    %Ub = [-Ub2;-Ub];
    %Lb = [-Lb2;-Lb];

end

function [dx] = odefun(x,u)
    b = 1.5 ; 
    L = 3 ;
    dx = [u(1)*cos(x(3))-(b/L)*u(1)*tan(u(2))*sin(x(3)) ;...
          u(1)*sin(x(3))+(b/L)*u(1)*tan(u(2))*cos(x(3)) ;...
          u(1)*tan(u(2))/L] ;
end

