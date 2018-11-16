data = load('./SimPkg_F18_V1/SimPkg_F18_V1/TestTrack.mat');
data = data.TestTrack;

plot(data.cline(1,:), data.cline(2,:), 'x')
hold on

%% ff gives a 3rd order polynomial  f = sum(a_i*x^(4-i)) i = 1 .. 4

for i = 1:length(data.cline(1,:))-1
    ff = csape(data.cline(1, i:i+1), data.cline(2,i:i+1), tan(data.theta(i:i+1)));
    if abs(ff.coefs(2)) > 1e-4 || abs(ff.coefs(1)) > 1e-4
        disp("not linear")
    end
    fnplt(ff, [data.cline(1,i), data.cline(1,i+1)])
    hold on
end
