%this function is designed as a cntroller for vehicle to perform feedback
%control alng the generated trjectory (with or whout obstacles)
%time step is 0.1s
%number of time steps between two trajectory points are dynamically
%determined
load TestTrack.mat
[Traj_ref_x, Traj_ref_y, Traj_ref_psi] = genTrajectory(TestTrack);
timestep = 0.1;