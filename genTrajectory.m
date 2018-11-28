% this is a trajectory generation script to directly generate reference
% trajectory from centerline of the road regardless of the resitriction of
% the left and right sides of the road
% for obstacles, this function use rrt to generate ramp to avoid collision
% traj(1:3,:) = x;y;psi
function traj = genTrajectory
%% gennerate trajectory without obstacles

