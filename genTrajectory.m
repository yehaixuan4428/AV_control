% this is a trajectory generation script to directly generate reference
% trajectory from centerline of the road regardless of the resitriction of
% the left and right sides of the road
% for obstacles, this function use rrt to generate ramp to avoid collision
% traj(1:3,:) = x;y;psi
function [Traj_ref_x, Traj_ref_y, Traj_ref_psi] = genTrajectory(cline,theta)
Obs_default = NaN;
if nargin == 1
    Obs = Obs_default;
end
%% gennerate trajectory without obstacles
Traj_ref_x = cline(1,:);
Traj_ref_y = cline(2,:);
Traj_ref_psi = theta;
end
