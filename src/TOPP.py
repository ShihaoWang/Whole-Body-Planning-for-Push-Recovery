"""Retime a path subject to kinematic constraints
=================================================

In this example, we will see how can we retime a generic spline-based
path subject to kinematic constraints. This is very simple to do with
`toppra`, as we shall see below. First import the library.

"""

import toppra as ta
import toppra.constraint as constraint
import toppra.algorithm as algo
import numpy as np
import matplotlib.pyplot as plt
import time

ta.setup_logging("INFO")

################################################################################
# We generate a path with some random waypoints.

def generate_new_problem(seed=9):
    # Parameters
    N_samples = 5
    dof = 7
    np.random.seed(seed)
    way_pts = np.random.randn(N_samples, dof)
    return (
        np.linspace(0, 1, 5),
        way_pts,
        10 + np.random.rand(dof) * 20,
        10 + np.random.rand(dof) * 2,
    )
ss, way_pts, vlims, alims = generate_new_problem()
import ipdb; ipdb.set_trace()
ss = np.linspace(0, 1, 6)
way_points = [[0.016217,-0.203080, 1.096496, 0.359483, -0.000354, 0.024674],
[-0.011923, -0.571062, 2.015125, -0.572965, -0.000772, 0.289342],
[-0.105746, -0.196983, 2.039504, -1.246179, 0.065274, -0.050439],
[-0.348617, -0.016003, 1.628022, -1.247505, 0.308215, -0.691632],
[-0.328162, 0.133297, 1.164121, -1.149354, 0.301942, -0.772231],
[-0.233429, 0.297723, 0.321921, -0.672444, 0.197631, -0.759877]]
way_pts = np.array(way_points)
vlims = 10 * np.ones(6)
alims = 100 * np.ones(6)

################################################################################
# Define the geometric path and two constraints.
path = ta.SplineInterpolator(ss, way_pts)
pc_vel = constraint.JointVelocityConstraint(vlims)
pc_acc = constraint.JointAccelerationConstraint(alims)

################################################################################
# We solve the parametrization problem using the
# `ParametrizeConstAccel` parametrizer. This parametrizer is the
# classical solution, guarantee constraint and boundary conditions
# satisfaction.

instance = algo.TOPPRA([pc_vel, pc_acc], path, parametrizer="ParametrizeConstAccel")
# instance = algo.TOPPRAsd([pc_vel, pc_acc], path)
# instance.set_desired_duration(0.6)

jnt_traj = instance.compute_trajectory(0.2, 0.5)

################################################################################
# The output trajectory is an instance of
# :class:`toppra.interpolator.AbstractGeometricPath`.
ts_sample = np.linspace(0, jnt_traj.duration, 100)
qs_sample = jnt_traj(ts_sample)
qds_sample = jnt_traj(ts_sample, 1)
qdds_sample = jnt_traj(ts_sample, 2)
fig, axs = plt.subplots(3, 1, sharex=True)
for i in range(path.dof):
    # plot the i-th joint trajectory
    axs[0].plot(ts_sample, qs_sample[:, i])
    axs[1].plot(ts_sample, qds_sample[:, i])
    axs[2].plot(ts_sample, qdds_sample[:, i])
axs[2].set_xlabel("Time (s)")
axs[0].set_ylabel("Position (rad)")
axs[1].set_ylabel("Velocity (rad/s)")
axs[2].set_ylabel("Acceleration (rad/s2)")
plt.show()


################################################################################
# Optionally, we can inspect the output.
instance.compute_feasible_sets()
instance.inspect()
