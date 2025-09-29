# fenics code for the thermal simulation of chip for the publication  for frequency is 3.5 GHZ
import sys
import meshio
from fenics import *
from dolfin import * 
from mshr import *
from numpy import loadtxt
from petsc4py import PETSc
import numpy as np
import time
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
T = 1200*3.125e-6                                                          # T is final time or total time (s)
num_steps = 1200                                                  # num_steps is the number of time steps
t = 0
Train_steps = num_steps
dt = T/num_steps                                                # time step size 
Rb = 2.26                                                        # the thermal resistor of convectional face
h_c = 2.40598e4
Ta = 0    
N_mode =12
ls = 200
ws = 200
hs = 17
mesh = BoxMesh(Point(0,0,0), Point(0.014,0.012,0.0002977976),ls-1,ws-1,hs-1)
coor1 = mesh.coordinates()
lr = coor1[:,0].max()
ll = coor1[:,0].min()
wb = coor1[:,1].min()
wt = coor1[:,1].max()
hmax = coor1[:,2].max()
hmin = coor1[:,2].min()
V = FunctionSpace(mesh, 'P', 1)
NU = 11
Num_nodes = mesh.num_vertices()
thick_actl =  0.0000557976
#number_mode = loadtxt('config_block.txt')
#define thermal conductivity 
#define thermal conductivity 
tol = 1E-14
k_0 = 100                                                                                                               # silicon conductivity      (w/(m.k))
k_1 = 1.2                                                   #oxide silicon thermal conductivity
kappa = Expression('x[2] <= 0.000242 + tol ? k_0 : k_1', degree=0,tol=tol, k_0=k_0, k_1=k_1) #define subdomain
#define density 
D0 = 2.33e3                                                         # silicon density    (kg/m^3)
D1 = 2.65e3                                                        # oxide silicon density
DS1 = Expression('x[2] <= 0.000242 + tol ? D0 : D1', degree=0,tol=tol, D0=D0, D1=D1) #define subdomain
#define specific heat
c0 = 751.1                                                         # silicon specific heat   (J/(kg.k))
c1 = 680                                                         # oxide silicon specific heat
sc = Expression('x[2] <= 0.000242 + tol ? c0 : c1', degree=0,tol=tol, c0=c0, c1=c1) #define subdomain

#define power source term 

T_integral = 0
#define initial value 
u0 = Constant(Ta)                                                                   # Ta is initial temperature 
u_n = interpolate(u0,V)
# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
solution = []
for n in range(0,Train_steps):
	solution.append('u1')
# Collect Neumann integrals
#a = DS1*sc*u*v*dx + dt*kappa*dot(grad(u), grad(v))*dx + sum(integrals_R_a)
coor = mesh.coordinates()
v2d = vertex_to_dof_map(V)
h = coor[:,2].max()
#CU = loadtxt('pod_result/CU.txt')
#print(len(CU[0]))
# if we have the solutiuon, we just need to  Load solution# #######################read the solution from a solution file ######################
u = Function(V)
n = 100
solution_load_file_name = "./solution_8_outof_11/file_" + str(n) + "h5"
solution_file = HDF5File(mesh.mpi_comm(), solution_load_file_name, "r")
solution_file.read(u, "solution")
solution[n] = interpolate(u0,V)
solution[n].assign(u)
solution_file.close()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~generate podmode~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ###########################################################################

CU = []
for n_t in range(0, 1):
    y = np.linspace(wb+ tol, wt- tol, 100)
    points = [(0.00361, y_,0.000242) for y_ in y]  # 2D points
    p_line = np.array([solution[n](point) for point in points])
    if n_t ==0:
        header_data = 'pod_result/y_block'
        data_file_name = header_data + '.txt'
        data_file = open(data_file_name,'w')
        for ii in range(0,len(y)):
            data_file.write('%.16g\n' % (y[ii]))
        data_file.close()
    header_data = 'pod_result/T_Y_'+str(n_t + 1)+'_mode'
    data_file_name = header_data + '.txt'
    data_file = open(data_file_name,'w')
    for ii in range(0,len(y)):
        data_file.write('%.16g\n' % (p_line[ii]))
    data_file.close()
    y = np.linspace(ll + tol, lr- tol, 100)
    points = [(y_,0.009205,0.000242) for y_ in y]
    p_line = np.array([solution[n](point) for point in points])
    if n_t == 0:
        header_data = 'pod_result/x_block'
        data_file_name = header_data + '.txt'
        data_file = open(data_file_name,'w')
        for ii in range(0,len(y)):
            data_file.write('%.16g\n' % (y[ii]))
        data_file.close()
    header_data = 'pod_result/T_X_'+str(n_t +1)+'_mode'
    data_file_name = header_data + '.txt'
    data_file = open(data_file_name,'w')
    for ii in range(0,len(y)):
        data_file.write('%.16g\n' % (p_line[ii]))
    data_file.close()

