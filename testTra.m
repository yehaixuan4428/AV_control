b = 1.45 ; 
L = 2.8;
dt=0.01;

global f_cl;
global f_bl;
global f_br;
global data;
[data, f_cl, f_bl, f_br] = genBounds('./SimPkg_F18_V1/SimPkg_F18_V1/TestTrack.mat');

nsteps = 50;
%global R;

%1.1
ub = [repmat([300; 100 ;pi], nsteps, 1); repmat([50;0.5], nsteps-1, 1)];
lb = [repmat([200; -200 ;-pi], nsteps, 1); repmat([0;-0.5], nsteps-1, 1)];

options = optimoptions('fmincon','SpecifyConstraintGradient',true,...
                        'FiniteDifferenceType', 'central', ...
                        'FiniteDifferenceStepSize', 1e-10, ...
                       'CheckGradients', true,...
                       'Display', 'iter');
                       %'MaxIter', 500);


cf=@costfun
nc=@nonlcon

x0=zeros(1,5*nsteps-2);
i = 1:nsteps;
xIndex = 3*i-2;
yIndex = xIndex+1;
psiIndex = xIndex+2;
%x0(1) = 287;
%x0(2) = -176.0;
%x0(3) = 2.0;
x0(xIndex) = linspace(287, data.cline(1, 4), nsteps);
x0(yIndex) = linspace(-176, data.cline(2, 4), nsteps);
x0(3) = 2.0;
x0(psiIndex(2:end)) = atan2(x0(yIndex(2:end)), x0(xIndex(2:end)));

%x0(xIndex(2:end)) = 287 + 50*cos(2.0)*(i(2:end) - 1);
%x0(yIndex(2:end)) = -176 + 50*sin(2.0)*(i(2:end) - 1);
%x0(psiIndex) = 2.0;

%figure
%plot(x0(xIndex), x0(yIndex), '-^')
%hold on
%plot(data.bl(1,1:15), data.bl(2,1:15), 'x', data.br(1,1:15), data.br(2,1:15), 'x')

i = 1:nsteps-1;
uIndex = nsteps*3 - 1 + 2*i;
deltaIndex = uIndex+1;
x0(uIndex) = 49.0;
x0(uIndex(1)) = 5.0;
Y = [];
U = [];

nTotalSeg = 1;
stepSet = [100, 100];

for nSeg = 1:nTotalSeg
    nSeg

    startState = [x0(1:3), x0(uIndex(1):deltaIndex(1))];
    z=fmincon(@(z) cf(z, [data.cline(:, 9);data.theta(9)], nsteps), x0,[],[],[],[],lb',ub',@(z) nc(z, nsteps, startState, func{nSeg}, nSeg),options);
    Y0=reshape(z(1:3*nsteps),3,nsteps)';
    U0=reshape(z(3*nsteps+1:end),2,nsteps-1);
    Y = [Y;Y0];
    U = [U,U0];
    %nsteps = stepSet(nSeg);
    i = 1:nsteps;
    xIndex = 3*i-2;
    yIndex = xIndex+1;
    psiIndex = xIndex+2;
    i = 1:nsteps-1;
    uIndex = nsteps*3 - 1 + 2*i;
    deltaIndex = uIndex+1;
    x0=zeros(1,5*nsteps-2);
    x0(1:3) = Y0(end,1:3);
    x0(uIndex(1):deltaIndex(1)) = U0(:, end)';
    ub = [repmat([1500; 900 ;pi], nsteps, 1); repmat([50;0.5], nsteps-1, 1)];
    lb = [repmat([200; -200 ;-pi], nsteps, 1); repmat([0;-0.5], nsteps-1, 1)];
end

%Y0=reshape(z(1:3*nsteps),3,nsteps)';
%U=reshape(z(3*nsteps+1:end),2,nsteps-1);
u=@(t) [interp1(0:dt:(nsteps-2)*dt,U(1,:),t,'previous','extrap');...
        interp1(0:dt:(nsteps-2)*dt,U(2,:),t,'previous','extrap')];
%Y1 = zeros(nsteps, 3);
%Y1(1,:) = Y0(1,:);
%for i = 1:nsteps - 1
    %Y1(i+1, :) = Y1(i,:) + dt.*(odefun(Y1(i,:), u(dt*(i-1))))';
%end
[T2,Y2]=ode45(@(t,x) odefun(x,u(t)),[0:dt:(nsteps-1)*dt], [287; -176.0; 2.0]);
%%plot(Y1(:,1),Y1(:,2));
figure
%plot(Y0(:,1),Y0(:,2),'x', Y1(:,1),Y1(:,2), 'o', Y2(:,1), Y2(:,2), '^');
%hold on
plot(Y(:, 1),Y(:,2), Y2(:,1), Y2(:,2), '--')
hold on
plot(data.bl(1,1:15), data.bl(2,1:15), 'x', data.br(1,1:15), data.br(2,1:15), 'x')

%%theta = 0:0.01:2*pi;
%%hold on
%%plot((0.7*cos(theta)+3.5),(0.7*sin(theta)-0.5))
%%hold on
%%plot(0,0,'x');
legend('fmincon', 'ode45 trajectory', 'left', 'right');
%%ylim([-2,2]);
%%xlim([-1,8]);
%xlabel('x');
%ylabel('y');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1.2
function [g,h,dg,dh]=nonlcon(z, nsteps, startStates, interpFunc, idx)
    b = 1.45; 
    L = 2.8;
    dt = 0.01;
    i = 1:nsteps-1;
    uIndex = nsteps*3 - 1 + 2*i;
    deltaIndex = uIndex+1;
    i = 1:nsteps;
    xIndex = 3*i-2;
    yIndex = xIndex+1;
    psiIndex = xIndex+2;
    dg = zeros(5*nsteps - 2, nsteps);

    if idx == 1
        g = interpFunc(z(xIndex))' - z(yIndex);
        for i = 1:nsteps
            dg(xIndex(i):yIndex(i), i) = [2*interpFunc.p1*xIndex(i) + interpFunc.p2; -1.0];
        end
    else
        g = interpFunc(z(yIndex))' - z(xIndex);
        for i = 1:nsteps
            dg(xIndex(i):yIndex(i), i) = [-1.0; 2*interpFunc.p1*yIndex(i) + interpFunc.p2];
        end
    end

    h = zeros(3*nsteps + 2,1);
    h(1:3) = z(xIndex(1):psiIndex(1));
    h(1:3) = h(1:3) - startStates(1:3)';

    h(xIndex(2:end)) = z(xIndex(2:end)) - z(xIndex(1:end-1)) - dt*(z(uIndex).*cos(z(psiIndex(1:end-1))) - b/L*z(uIndex).*tan(z(deltaIndex)).*sin(z(psiIndex(1:end-1))));
    h(yIndex(2:end)) = z(yIndex(2:end)) - z(yIndex(1:end-1)) - dt*(z(uIndex).*sin(z(psiIndex(1:end-1))) + b/L*z(uIndex).*tan(z(deltaIndex)).*cos(z(psiIndex(1:end-1))));
    h(psiIndex(2:end)) = z(psiIndex(2:end)) - z(psiIndex(1:end-1)) - dt/L*z(uIndex).*tan(z(deltaIndex));

    h(psiIndex(end) + 1) = z(uIndex(1)) - startStates(4);
    h(psiIndex(end) + 2) = z(deltaIndex(1)) - startStates(5);

    dh = zeros(5*nsteps - 2, 3*nsteps+2);
    for i = 1:nsteps -1
        dh(1:3,1:3) = eye(3);
        dh(xIndex(i):psiIndex(i+1), xIndex(i+1):psiIndex(i+1)) = [ ...
    -1.0, 0.0, 0.0;
    0.0, -1.0, 0.0;
    -dt*(-z(uIndex(i))*sin(z(psiIndex(i)))-b/L*z(uIndex(i))*tan(z(deltaIndex(i)))*cos(z(psiIndex(i)))), -dt*(z(uIndex(i))*cos(z(psiIndex(i)))-b/L*z(uIndex(i))*tan(z(deltaIndex(i)))*sin(z(psiIndex(i)))), -1.0;
    eye(3) ];
        dh(uIndex(i):deltaIndex(i), xIndex(i+1):psiIndex(i+1)) = -dt*[
                cos(z(psiIndex(i))) - b/L*tan(z(deltaIndex(i)))*sin(z(psiIndex(i))),  sin(z(psiIndex(i))) + b/L*tan(z(deltaIndex(i)))*cos(z(psiIndex(i))), 1.0/L*tan(z(deltaIndex(i)));
                -b/L*z(uIndex(i))/(cos(z(deltaIndex(i))))^2*sin(z(psiIndex(i))), b/L*z(uIndex(i))/(cos(z(deltaIndex(i))))^2*cos(z(psiIndex(i))), z(uIndex(i))/L/(cos(z(deltaIndex(i))))^2];
    end
    dh(uIndex(1), 3*nsteps+1) = 1.0;
    dh(deltaIndex(1), end) = 1.0;

    %%% first term
    %dh1 = eye(603, 363);
    %%% second term first component
    %value = -ones(1, 120*3);
    %dh2 = sparse([xIndex(1:end-1) yIndex(1:end-1) psiIndex(1:end-1)], [xIndex(2:end) yIndex(2:end) psiIndex(2:end)], value, 603,363);

    %%% third term first component
    %row = repmat([psiIndex(1:end-1) uIndex deltaIndex], 1, 2);
    %row = [row uIndex deltaIndex];
    %col = [repmat(xIndex(2:end), 1, 3) repmat(yIndex(2:end), 1, 3) repmat(psiIndex(2:end), 1, 2)];
    %value = -dt*[
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

function [J, dJ] = costfun(z, endPoint, nsteps)
    xIndex = 1:3:nsteps*3;
    yIndex = 2:3:nsteps*3;
    psiIndex = 3:3:nsteps*3;
    uIndex = nsteps*3+1:2:nsteps*5-2;
    deltaIndex = nsteps*3+2:2:nsteps*5-2;
    %J = sum(  1e-5*((z(xIndex) - endPoint(1)).^2 + (z(yIndex) - endPoint(2)).^2) + (z(psiIndex) - endPoint(3)).^2) + sum(z(deltaIndex).^2);
    J = sum(  1e-5*((z(xIndex) - endPoint(1)).^2 + (z(yIndex) - endPoint(2)).^2) + (z(psiIndex) - endPoint(3)).^2);
    dJ = 2.0*z';
    dJ(xIndex) = 1e-5*dJ(xIndex) - 1e-5*endPoint(1)*2;
    dJ(yIndex) = 1e-5*dJ(yIndex) - 1e-5*endPoint(2)*2;
    dJ(psiIndex) = dJ(psiIndex) - endPoint(3)*2;
    dJ(uIndex) = 0.0;
    dJ(deltaIndex) = 0.0;
end


function [dx] = odefun(x,u)
    b = 1.45 ; 
    L = 2.8 ;
    dx = [u(1)*cos(x(3))-(b/L)*u(1)*tan(u(2))*sin(x(3)) ;...
          u(1)*sin(x(3))+(b/L)*u(1)*tan(u(2))*cos(x(3)) ;...
          u(1)*tan(u(2))/L];
end
