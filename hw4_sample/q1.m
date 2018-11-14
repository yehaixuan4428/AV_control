b = 1.5 ; 
L = 3 ;
dt=0.05;
nsteps = 121;

%remember the format for z is as follows:
%z=[x0 y0 th0 x1 y1 th1 ... xn yn thn u0 d0 ... u(n-1) d(n-1)]';
    
%1.1
ub = [repmat([8;3;pi/2.0], 121, 1); repmat([1;0.5], 120, 1)];
lb = [repmat([-1;-3;-pi/2.0], 121, 1); repmat([0;-0.5], 120, 1)];

%1.4
%%%%%%%%%%%%%%%% no need to change these lines  %%%%%%%%%%%%%%%%%%%%%%
options = optimoptions('fmincon','SpecifyConstraintGradient',true,...
                       'SpecifyObjectiveGradient',true, 'Display', 'iter');
x0=zeros(1,5*nsteps-2);
cf=@costfun
nc=@nonlcon
z=fmincon(cf,x0,[],[],[],[],lb',ub',nc,options);

Y0=reshape(z(1:3*nsteps),3,nsteps)';
U=reshape(z(3*nsteps+1:end),2,nsteps-1);
u=@(t) [interp1(0:dt:119*dt,U(1,:),t,'previous','extrap');...
        interp1(0:dt:119*dt,U(2,:),t,'previous','extrap')];
[T1,Y1]=ode45(@(t,x) odefun(x,u(t)),[0:dt:120*dt],[0 0 0]);
[T2,Y2]=ode45(@(t,x) odefun(x,u(t)),[0:dt:120*dt],[0 0 -0.01]);
plot(Y0(:,1),Y0(:,2),Y1(:,1),Y1(:,2),Y2(:,1),Y2(:,2))
theta = 0:0.01:2*pi;
hold on
plot((0.7*cos(theta)+3.5),(0.7*sin(theta)-0.5))
hold on
plot(0,0,'x');
legend('fmincon trajectory','ode45 trajectory using x0 = [0;0;0]',...
    'ode45 trajectory using x0 = [0;0;-0.01]','Buffered Obstacle','Start');
ylim([-2,2]);
xlim([-1,8]);
xlabel('x');
ylabel('y');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1.2
function [g,h,dg,dh]=nonlcon(z)
    b = 1.5;
    L = 3;
    i = 1:120;
    uIndex = 362 + 2*i;
    deltaIndex = uIndex+1;
    i = 1:121;
    xIndex = 3*i-2;
    yIndex = xIndex+1;
    psiIndex = xIndex+2;

    g = 0.7^2 - (z(xIndex) - 3.5).^2 - (z(yIndex) + 0.5).^2;
    dg = zeros(603, 121);

    h = zeros(363,1);
    h(xIndex(2:end)) = z(xIndex(2:end)) - z(xIndex(1:end-1)) - 0.05*(z(uIndex).*cos(z(psiIndex(1:end-1))) - b/L*z(uIndex).*tan(z(deltaIndex)).*sin(z(psiIndex(1:end-1))));
    h(yIndex(2:end)) = z(yIndex(2:end)) - z(yIndex(1:end-1)) - 0.05*(z(uIndex).*sin(z(psiIndex(1:end-1))) + b/L*z(uIndex).*tan(z(deltaIndex)).*cos(z(psiIndex(1:end-1))));
    h(psiIndex(2:end)) = z(psiIndex(2:end)) - z(psiIndex(1:end-1)) - 0.05/L*z(uIndex).*tan(z(deltaIndex));

    dh = zeros(603, 363);

    for i = 1:120
        dg(xIndex(i), i) = -2.0*(z(xIndex(i)) - 3.5);
        dg(yIndex(i), i) = -2.0*(z(yIndex(i)) + 0.5);
        dh(1:3,1:3) = eye(3);
        dh(xIndex(i):psiIndex(i+1), xIndex(i+1):psiIndex(i+1)) = [ ...
    -1.0, 0.0, 0.0;
    0.0, -1.0, 0.0;
    -0.05*(-z(uIndex(i))*sin(z(psiIndex(i)))-b/L*z(uIndex(i))*tan(z(deltaIndex(i)))*cos(z(psiIndex(i)))), -0.05*(z(uIndex(i))*cos(z(psiIndex(i)))-b/L*z(uIndex(i))*tan(z(deltaIndex(i)))*sin(z(psiIndex(i)))), -1.0;
    eye(3) ];
        dh(uIndex(i):deltaIndex(i), xIndex(i+1):psiIndex(i+1)) = -0.05*[
                cos(z(psiIndex(i))) - b/L*tan(z(deltaIndex(i)))*sin(z(psiIndex(i))),  sin(z(psiIndex(i))) + b/L*tan(z(deltaIndex(i)))*cos(z(psiIndex(i))), 1.0/L*tan(z(deltaIndex(i)));
                -b/L*z(uIndex(i))/(cos(z(deltaIndex(i))))^2*sin(z(psiIndex(i))), b/L*z(uIndex(i))/(cos(z(deltaIndex(i))))^2*cos(z(psiIndex(i))), z(uIndex(i))/L/(cos(z(deltaIndex(i))))^2];
    end
    dg(xIndex(end), end) = -2.0*(z(xIndex(end)) - 3.5);
    dg(yIndex(end), end) = -2.0*(z(yIndex(end)) + 0.5);

    %%% first term
    %dh1 = eye(603, 363);
    %%% second term first component
    %value = -ones(1, 120*3);
    %dh2 = sparse([xIndex(1:end-1) yIndex(1:end-1) psiIndex(1:end-1)], [xIndex(2:end) yIndex(2:end) psiIndex(2:end)], value, 603,363);

    %%% third term first component
    %row = repmat([psiIndex(1:end-1) uIndex deltaIndex], 1, 2);
    %row = [row uIndex deltaIndex];
    %col = [repmat(xIndex(2:end), 1, 3) repmat(yIndex(2:end), 1, 3) repmat(psiIndex(2:end), 1, 2)];
    %value = -0.05*[
                %-z(uIndex).*sin(z(psiIndex(1:end-1))) - b/L*z(uIndex).*tan(z(deltaIndex)).*cos(z(psiIndex(1:end-1))) ...
                %cos(z(psiIndex(1:end-1))) - b/L*tan(z(deltaIndex)).*sin(z(psiIndex(1:end-1))) ...
                %-b/L*z(uIndex)./(cos(z(deltaIndex))).^2.*sin(psiIndex(1:end-1)) ...
                %z(uIndex).*cos(psiIndex(1:end-1)) - b/L*z(uIndex).*tan(z(deltaIndex)).*sin(z(psiIndex(1:end-1))) ...
                %sin(z(psiIndex(1:end-1))) + b/L*tan(z(deltaIndex)).*cos(z(psiIndex(1:end-1))) ...
                %b/L*z(uIndex)./(cos(z(deltaIndex))).^2.*cos(psiIndex(1:end-1)) ...
                %1.0/L*tan(z(deltaIndex)) ...
                %z(uIndex)./(cos(z(deltaIndex))).^2 ...
            %];
    %dh3 = sparse(row, col, value, 603, 363);
    %dh = dh1 + dh2 + dh3;

    % size of g must be 121 x 1 (no.of time steps);
    % size of dg must be 603 x 121 = Transpose(no. of time steps x no. of elements in 'z');
    % size of h must be 363  1 ((no. of time steps * no. of states) x 1)
    % size of dh must be 603 x 363 = Transpose((no. of time steps * no. of states) x no. of elements in 'z') ;
end

%1.3q
function [J, dJ] = costfun(z)
    % size of J must be 1 x 1
    % size of dJ must be 1 x 603 (1 x no. of elements in 'z')
    xIndex = 1:3:363;
    yIndex = 2:3:363;
    psiIndex = 3:3:363;
    uIndex = 364:2:603;
    deltaIndex = 365:2:603;

    J = sum((z(xIndex) - 7.0).^2 + z(yIndex).^2 + z(psiIndex).^2) + sum(z(uIndex).^2 + z(deltaIndex).^2);
    dJ = 2.0*z';
    dJ(xIndex) = dJ(xIndex) - 14.0;

end

function [dx] = odefun(x,u)
    b = 1.5 ; 
    L = 3 ;
    dx = [u(1)*cos(x(3))-(b/L)*u(1)*tan(u(2))*sin(x(3)) ;...
          u(1)*sin(x(3))+(b/L)*u(1)*tan(u(2))*cos(x(3)) ;...
          u(1)*tan(u(2))/L] ;
end
