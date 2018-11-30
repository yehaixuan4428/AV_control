% this is a trajectory generation script to directly generate reference
% trajectory from centerline of the road regardless of the resitriction of
% the left and right sides of the road
% for obstacles, this function use rrt to generate ramp to avoid collision
% traj(1:3,:) = x;y;psi
function traj = genTrajectory(TestTrack,Obs)
Obs_default = NaN;
if nargin == 1
    Obs = Obs_default;
end
%% gennerate trajectory without obstacles
traj(1,:) = TestTrack.cline(1,:);
traj(2,:) = TestTrack.cline(2,:);
traj(3,:) = TestTrack.theta;
end
