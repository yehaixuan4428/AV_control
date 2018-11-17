function [a,b,c]=linearInterpolate(T, index, cond)
%cond=1 gives left bound, cond=2 gives right bound
%T: one of the bound of TestTrack.mat data struct
%index: index of the closest bound's index of the TestTrack

if cond == 1
    t=T.bl;
elseif cond == 2
    t=T.br;
end

%find by using csape
%!!!The 3rd input of csape do not effect anything
f = csape(t(1, index:index+1), t(2,index:index+1), tan(T.theta(index:index+1)));
f_test = csape(t(1, index:index+1), t(2,index:index+1));

a = f.coefs(3);
b = f.coefs(4);
c = a*t(1,index)+b*t(2,index);

a_test = f_test.coefs(3);
b_test = f_test.coefs(4);
a==a_test
b==b_test

%find by linear interpolate
A = [t(1,index:index+1)',t(2,index:index+1)'];
sol=linsolve(A,[1;1]);
a_lin=sol(1);
b_lin=sol(2);
c_lin=1;

figure(1)
fnplt(f, [t(1, index), t(1, index+1)]);
figure(2)
x = min(t(1,index),t(1,index+1)):0.01:max(t(1,index),t(1,index+1));
y = (c_lin-a_lin*x)./b_lin;
plot(x,y);
end