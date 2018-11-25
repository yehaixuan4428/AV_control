function [a,b,c,func] = CurveFitting(bd)
%a,b,c are column vectors of the coefficients of f(x)=ax^2+bx+c 
%func is a cell array cointains all fitted functions
%{
Points 10~15 are unable to be directly fitted by a 2nd 
order polynomial.
In order to get a non-linear differentiable fitted functions,
their x and y values are swaped.
%}
%bd = TestTrack.lb(rb) 
f=fit(bd(1,1:9)',bd(2,1:9)','poly2');
a = [f.p1];
b = [f.p2];
c = [f.p3];
func = {f};
f=fit(bd(2,10:15)',bd(1,10:15)','poly2');
a = [a;f.p1];
b = [b;f.p2];
c = [c;f.p3];
func = [func;{f}];
interval = 6; %fit every k points  
i=16;
while (i+interval-1)<246 
    f=fit(bd(1,i:i+interval-1)',bd(2,i:i+interval-1)','poly2');
    a = [a;f.p1];
    b = [b;f.p2];
    c = [c;f.p3];
    func = [func;{f}];
    i=i+interval;
end
f=fit(bd(1,i:end)',bd(1,i:end)','poly2');
a = [a;f.p1];
b = [b;f.p2];
c = [c;f.p3];
func = [func;{f}];
end