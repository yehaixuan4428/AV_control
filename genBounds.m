function [data, f_cl, f_bl, f_br] = genBounds(fileName)
    data = load(fileName);
    data = data.TestTrack;
    nSeg = length(data.cline(1,:)) - 1;
    %% ff gives a 3rd order polynomial  f = sum(a_i*x^(4-i)) i = 1 .. 4

    for i = 1:nSeg
        f = csape(data.cline(1, i:i+1), data.cline(2,i:i+1), tan(data.theta(i:i+1)));
        f_cl(i,:) = f.coefs(3:4);
        f = csape(data.bl(1, i:i+1), data.bl(2,i:i+1), tan(data.theta(i:i+1)));
        f_bl(i,:) = f.coefs(3:4);
        f = csape(data.br(1, i:i+1), data.br(2,i:i+1), tan(data.theta(i:i+1)));
        f_br(i,:) = f.coefs(3:4);

    %%% Check if the polynomial is linear
    %if abs(f_cl(i).coefs(2)) > 1e-4 || abs(f_cl(i).coefs(1)) > 1e-4
        %disp("not linear")
    %end
    %if abs(f_bl(i).coefs(2)) > 1e-4 || abs(f_bl(i).coefs(1)) > 1e-4
        %disp("not linear")
    %end
    %if abs(f_br(i).coefs(2)) > 1e-4 || abs(f_br(i).coefs(1)) > 1e-4
        %disp("not linear")
    %end
    end
end


%%% plot
%for i = 1:nSeg
    %fnplt(f_cl(i), [data.cline(1, i), data.cline(1, i+1)]);
    %hold on
    %fnplt(f_bl(i), [data.bl(1, i), data.bl(1, i+1)]);
    %hold on
    %fnplt(f_br(i), [data.br(1, i), data.br(1, i+1)]);
    %hold on
%end
%plot(data.cline(1,:), data.cline(2,:), 'x')
%hold on
%plot(data.bl(1,:), data.bl(2,:), 'x')
%hold on
%plot(data.br(1,:), data.br(2,:), 'x')
%hold on

