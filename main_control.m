%this function is designed as a cntroller for vehicle to perform feedback
%control alng the generated trjectory (with or whout obstacles)
%time step is 0.1s
%number of time steps between two trajectory points are dynamically
%determined
clear
clc
load TestTrack.mat
Pcontroller = [0 0.98; 300 -5];
[Traj_ref_x, Traj_ref_y, Traj_ref_psi] = genTrajectory(TestTrack.cline,TestTrack.theta);
timestep = 0.1;
plot(TestTrack.bl(1,:),TestTrack.bl(2,:),TestTrack.br(1,:),TestTrack.br(2,:),Traj_ref_x,Traj_ref_y,'o','MarkerSize',2)
hold on
% define initial state
z(:,1) = [287;5;-176;0;2;0];
unom = 5;
u(:,1) = Pcontroller*([unom;Traj_ref_psi(2)]-...
        [z(2,1);z(5,1)]);
for i = 2:length(Traj_ref_x)
    while ~pt_threshold([z(1,end);z(3,end)],[Traj_ref_x(i);Traj_ref_y(i)],[Traj_ref_x(i-1);Traj_ref_y(i-1)])
        Pcontroller = [0 0.98; 300 -5];
        if u(1,end)>0.5
            u(1,end) = 0.5;
        elseif u(1,end)<-0.5
            u(1,end) = -0.5;
        end
        if u(2,end)>2500
            u(2,end) = 2500;
        elseif u(2,end)<-5000
            u(2,end) = -5000;
        end
        dz = vehicle_model(z(:,end),u(:,end));
        z(:,size(z,2)+1) = z(:,end)+timestep*dz;
        a(size(z,2)+1) = Traj_ref_psi(i)-z(5,end);
        if Traj_ref_psi(i)-z(5,end)>0.18
            Pcontroller(1,2) = 20;
        end
        u(:,size(u,2)+1) = Pcontroller*([unom;Traj_ref_psi(i)]-...
        [z(2,end);z(5,end)]);
        plot(z(1,end),z(3,end),'.')
        hold on
        axis([200 1600 -200 1000])
    end
end

function k = pt_threshold(x,x_d,x_d_b)
dist = norm(x_d-x);
dist_d = 0.3*norm(x_d-x_d_b);
    if dist< dist_d
        k = true;
    else
        k = false;
    end
end