%this function is designed as a cntroller for vehicle to perform feedback
%control alng the generated trjectory (with or whout obstacles)
%time step is 0.1s
%number of time steps between two trajectory points are dynamically
%determined
clear
clc
load TestTrack.mat
load Pcontroller.mat
[Traj_ref_x, Traj_ref_y, Traj_ref_psi] = genTrajectory(TestTrack.cline,TestTrack.theta);
timestep = 0.1;
plot(Traj_ref_x,Traj_ref_y)
hold on
% define initial state
z(:,1) = [287;5;-176;0;2;0];
u(:,1) = Pcontroller*([Traj_ref_x(2);Traj_ref_y(2);Traj_ref_psi(2)]-...
        [z(1,1);z(3,1);z(5,1)]);
for i = 2:length(Traj_ref_x)
    while ~pt_threshold([z(1,end);z(3,end)],[Traj_ref_x(i);Traj_ref_y(i)],[Traj_ref_x(i-1);Traj_ref_y(i-1)])
        if u(2,end)>0.5
            u(2,end) = 0.5;
        elseif u(2,end)<-0.5
            u(2,end) = -0.5;
        end
        dz = vehicle_model(z(:,end),u(:,end));
        z(:,size(z,2)+1) = z(:,end)+timestep*dz;
        u(:,size(u,2)+1) = Pcontroller*([Traj_ref_x(i);Traj_ref_y(i);Traj_ref_psi(i)]-...
        [z(1,end);z(3,end);z(5,end)]);
        plot(z(1,end),z(3,end),'*')
        hold on
    end
end

function k = pt_threshold(x,x_d,x_d_b)
dist = norm(x_d-x);
dist_d = 0.025*norm(x_d-x_d_b);
    if dist< dist_d
        k = true;
    else
        k = false;
    end
end