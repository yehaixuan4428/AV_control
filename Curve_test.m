clear all; clc;
data = load('./SimPkg_F18_V1/SimPkg_F18_V1/TestTrack.mat');
data = data.TestTrack;
[a,b,c,func]=CurveFitting(data.bl);

%generate random points
x = 200+1300*rand([1,200]);
y = -200+1000*rand([1,200]);
seg = length(a);

figure() %Test the 1~9 points
hold on
f = func(1);
f = f{1};
% f=fit(data.bl(1,1:9)',data.bl(2,1:9)','poly2');
plot(f,data.bl(1,:)',data.bl(2,:)')
axis([200 1500 -200 800])
for i=1:200
    if y(i)-f(x(i))>0
        scatter(x(i),y(i),'g')
    else
        scatter(x(i),y(i),'r')
    end
end

figure() %Test the 10~15 points
plot(data.bl(1,:)',data.bl(2,:)')
axis([200 1500 -200 800])
hold on
f = func(2);
f = f{1};
for i=1:200
    if x(i)-f(y(i))>0
        scatter(x(i),y(i),'g')
    else
        scatter(x(i),y(i),'r')
    end
end


figure() %Test the 16~21 points
hold on
f = func(3);
f = f{1};
% f=fit(data.bl(1,16:21)',data.bl(2,16:21)','poly2');
plot(f,data.bl(1,:)',data.bl(2,:)')
axis([200 1500 -200 800])

for i=1:200
    if y(i)-f(x(i))<0
        scatter(x(i),y(i),'g')
    else
        scatter(x(i),y(i),'r')
    end
end


figure() %Test the 22~27 points
hold on
f = func(4);
f = f{1};
% f=fit(data.bl(1,22:27)',data.bl(2,22:27)','poly2');
plot(f,data.bl(1,:)',data.bl(2,:)')
axis([200 1500 -200 800])

for i=1:200
    if y(i)-f(x(i))<0
        scatter(x(i),y(i),'g')
    else
        scatter(x(i),y(i),'r')
    end
end
