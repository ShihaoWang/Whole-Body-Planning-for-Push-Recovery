import sys, os, time
from klampt import *
from klampt import vis
from klampt.vis.glrobotprogram import GLSimulationPlugin
from klampt.model.trajectory import Trajectory, RobotTrajectory
import ipdb
import copy
from scipy.spatial import ConvexHull
import draw_hull
from OpenGL.GL import *
import math
import numpy as np
import random

class MyGLPlugin(vis.GLPluginInterface):
    def __init__(self, world):
        vis.GLPluginInterface.__init__(self)
        self.world = world
        self.quit = False
        self.starp = False

    def mousefunc(self, button, state, x, y):
        print("mouse",button,state,x,y)
        if button==2:
            if state==0:
                print("Click list...",[o.getName() for o in self.click_world(x,y)])
            return True
        return False

    def motionfunc(self, x, y, dx, dy):
        return False

    def keyboardfunc(self, c, x, y):
        print("Pressed", c)
        return True

    def click_world(self, x, y):
        """Helper: returns a list of world objects sorted in order of
        increasing distance."""
        #get the viewport ray
        (s, d) = self.click_ray(x, y)

        #run the collision tests
        collided = []
        for g in self.collider.geomList:
            (hit, pt) = g[1].rayCast(s, d)
            if hit:
                dist = vectorops.dot(vectorops.sub(pt, s), d)
                collided.append((dist,g[0]))
        return [g[1] for g in sorted(collided)]

def ContactDataUnplot(vis, ReachableContacts_data):
    RowNo, ColumnNo = ReachableContacts_data.shape
    RowStart = 0
    RowEnd = RowNo
    for i in range(RowStart, RowEnd):
        vis.remove("Point:" + str(i))

def ContactDataPlot(vis, ReachableContacts_data):
    RowNo, ColumnNo = ReachableContacts_data.shape
    RowStart = 0
    RowEnd = RowNo

    for i in range(RowStart, RowEnd):
        point_start = [0.0, 0.0, 0.0]
        ReachableContact_i = ReachableContacts_data[i]
        point_start[0] = ReachableContact_i[0]
        point_start[1] = ReachableContact_i[1]
        point_start[2] = ReachableContact_i[2]

        vis.add("Point:" + str(i), point_start)
        vis.hideLabel("Point:" + str(i), True)
        vis.setColor("Point:" + str(i),65.0/255.0, 199.0/255.0, 244.0/255.0, 1.0)

def WeightedContactDataPlot(vis, OptimalContact_data, OptimalContactWeights_data):
    scale = 1.0
    for i in range(OptimalContact_data.size/3):
        point_start = [0.0, 0.0, 0.0]
        ReachableContact_i = OptimalContact_data[i]
        point_start[0] = ReachableContact_i[0]
        point_start[1] = ReachableContact_i[1]
        point_start[2] = ReachableContact_i[2]

        point_end = [0.0, 0.0, 0.0]
        ReachableContactWeight_i = OptimalContactWeights_data[i]
        point_end[0] = point_start[0] + scale * ReachableContactWeight_i[0]
        point_end[1] = point_start[1] + scale * ReachableContactWeight_i[1]
        point_end[2] = point_start[2] + scale * ReachableContactWeight_i[2]

        vis.add("PointWeights:" + str(i), Trajectory([0, 1], [point_start, point_end]))
        vis.hideLabel("PointWeights:" + str(i), True)
        vis.setColor("PointWeights:" + str(i), 0.0, 204.0/255.0, 0.0, 1.0)
        vis.setAttribute("PointWeights:" + str(i), 'width', 5.0)

def WeightedContactDataUnPlot(vis, OptimalContact_data):
    for i in range(OptimalContact_data.size/3):
        vis.remove("PointWeights:" + str(i))

def Robot_Config_Plot(world, DOF, config_init):
    robot_viewer = MyGLPlugin(world)
    vis.pushPlugin(robot_viewer)
    vis.add("world", world)
    vis.show()

   # Here we would like to read point cloud for visualization of planning.
    # 1. All Reachable Points
    # IdealReachableContacts_data = ContactDataLoader("IdealReachableContact")
    # 2. Active Reachable Points
    # ActiveReachableContacts_data = ContactDataLoader("ActiveReachableContact")
    # 3. Contact Free Points
    ContactFreeContacts_data = ContactDataLoader("ContactFreeContact")
    # 4. Supportive Points
    SupportContacts_data = ContactDataLoader("SupportContact")
    # 5. Optimal Point
    OptimalContact_data = ContactDataLoader("OptimalContact")

    OptimalContactWeights_data = ContactDataLoader("OptimalContactWeights")
    # 6.
    TransitionPoints_data = ContactDataLoader("TransitionPoints")
    # 7.
    InitialTransitionPoints_data = ContactDataLoader("InitialTransitionPoints")
    # 8.
    ShiftedTransitionPoints_data = ContactDataLoader("ShiftedTransitionPoints")

    ReducedOptimalContact_data = ContactDataLoader("ReducedOptimalContact")

    ContactChoice = ShiftedTransitionPoints_data
    SimRobot = world.robot(0)
    SimRobot.setConfig(config_init)
    # import ipdb; ipdb.set_trace()
    while vis.shown():
        # This is the main plot program
        vis.lock()
        SimRobot.setConfig(config_init)
        WeightedContactDataPlot(vis, OptimalContact_data, OptimalContactWeights_data)
        ContactDataPlot(vis, ContactChoice)
        vis.unlock()
        time.sleep(0.1)
        WeightedContactDataUnPlot(vis, OptimalContact_data)
        ContactDataUnplot(vis, ContactChoice)

def RobotCOMPlot(SimRobot, vis):
    COMPos_start = SimRobot.getCom()
    COMPos_end = COMPos_start[:]
    COMPos_end[2] = COMPos_end[2] - 7.50
    vis.add("COM", Trajectory([0, 1], [COMPos_start, COMPos_end]))
    vis.hideLabel("COM",True)
    vis.setColor("COM", 0.0, 204.0/255.0, 0.0, 1.0)
    vis.setAttribute("COM",'width', 5.0)

def ContactDataGene(ReachableContacts_data):
    RowNo, ColumnNo = ReachableContacts_data.shape
    RowStart = 0
    RowEnd = RowNo

    point = []
    for i in range(RowStart, RowEnd):
        point_start = [0.0, 0.0, 0.0]
        ReachableContact_i = ReachableContacts_data[i]
        point_start[0] = ReachableContact_i[0]
        point_start[1] = ReachableContact_i[1]
        point_start[2] = ReachableContact_i[2]
        point.append(point_start)
    return ContactDataGene



def Configuration_Loader_fn(Config_Name):
    # This function is only used to load in the initial configuraiton
    # The initial file will be in the .config format
    with open(Config_Name,'r') as robot_angle_file:
        robotstate_angle_i = robot_angle_file.readlines()
    config_temp = [x.replace('\t',' ') for x in robotstate_angle_i]
    config_temp = [x.replace('\n','') for x in config_temp]
    config_temp = [float(i) for i in config_temp[0].split()]

    DOF = int(config_temp[0])
    # Config_Init = np.array(config_temp[1:])
    Config_Init = config_temp[1:]
    return DOF, Config_Init

def ContactDataLoader(IdealReachableContact):
    IdealReachableContacts = "/home/motion/Desktop/Whole-Body-Planning-for-Push-Recovery/build/" + IdealReachableContact + ".bin"
    f_IdealReachableContacts = open(IdealReachableContacts, 'rb')
    IdealReachableContacts_data = np.fromfile(f_IdealReachableContacts, dtype=np.double)
    IdealReachableContacts_data = IdealReachableContacts_data.reshape((IdealReachableContacts_data.size/3, 3))
    return IdealReachableContacts_data

def main(*arg):
    import ipdb; ipdb.set_trace()
    Robot_Option = "../user/"
    world = WorldModel()                    	# WorldModel is a pre-defined class
    # The next step is to load in robot's XML file
    XML_path = "/home/motion/Desktop/Whole-Body-Planning-for-Push-Recovery-Data/result/flat_1Contact/Environment.xml"
    result = world.readFile(XML_path)         	# Here result is a boolean variable indicating the result of this loading operation
    if not result:
        raise RuntimeError("Unable to load model " + XML_path)
    # In this case, what we have is a config
    ConfigName = "/home/motion/Desktop/Whole-Body-Planning-for-Push-Recovery-Data/result/flat_1Contact/26/InitConfig.config"
    DOF, Config_Init = Configuration_Loader_fn(ConfigName)
    Robot_Config_Plot(world, DOF, Config_Init)
if __name__ == "__main__":
    main("Config")
