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
Index_power = int(sys.argv[1])
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
T = 4800*3.125e-6                                                          # T is final time or total time (s)
num_steps = 4800                                                  # num_steps is the number of time steps
t = 0
Train_steps = num_steps
dt = T/num_steps                                                # time step size 
Rb = 2.26                                                        # the thermal resistor of convectional face
h_c = 2.40598e4
Ta = 0    
N_mode =12
ls = 150
ws = 150
hs = 17
mesh = BoxMesh(Point(0,0,0), Point(0.014,0.012,7.2e-4),ls-1,ws-1,hs-1)
coor1 = mesh.coordinates()
lr = coor1[:,0].max()
ll = coor1[:,0].min()
wb = coor1[:,1].min()
wt = coor1[:,1].max()
hmax = coor1[:,2].max()
hmin = coor1[:,2].min()
V = FunctionSpace(mesh, 'P', 1)
NU = 13
Num_nodes = mesh.num_vertices()
thick_actl =  1.35e-4
#number_mode = loadtxt('config_block.txt')
#define thermal conductivity 
#define thermal conductivity 
tol = 1E-14
k_0 = 100                                                                                                               # silicon conductivity      (w/(m.k))
k_1 = 1.2                                                   #oxide silicon thermal conductivity
kappa = Expression('x[2] <= 5.85e-4 + tol ? k_0 : k_1', degree=0,tol=tol, k_0=k_0, k_1=k_1) #define subdomain
#define density 
D0 = 2.33e3                                                         # silicon density    (kg/m^3)
D1 = 2.65e3                                                        # oxide silicon density
DS1 = Expression('x[2] <= 5.85e-4 + tol ? D0 : D1', degree=0,tol=tol, D0=D0, D1=D1) #define subdomain
#define specific heat
c0 = 751.1                                                         # silicon specific heat   (J/(kg.k))
c1 = 680                                                         # oxide silicon specific heat
sc = Expression('x[2] <= 5.85e-4 + tol ? c0 : c1', degree=0,tol=tol, c0=c0, c1=c1) #define subdomain

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
n = Train_steps -1
solution_load_file_name = "./solution_13blocks_pre"+str(Index_power)+"/file_" + str(n) + "h5"
solution_file = HDF5File(mesh.mpi_comm(), solution_load_file_name, "r")
solution_file.read(u, "solution")
solution[n] = interpolate(u0,V)
solution[n].assign(u)
solution_file.close()
"""
vtkfile_name =  'FEM_pre'+str(Index_power)+'.pvd'
vtkfile = File(vtkfile_name)
vtkfile << u
"""

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~generate podmode~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ###########################################################################
obj = {}
for i in [2]:
#for i in range(0,NU):
    obj['podmode' + str(i)] = []
    for n in range(0,N_mode):
        obj['podmode' + str(i)].append('u1')
        podmode_load_file_name = "./podmode_" + str(i)+"/mode_" + str(n) + "h5"
        podmode_file = HDF5File(mesh.mpi_comm(), podmode_load_file_name, "r")
        podmode_file.read(u, "solution")
        obj['podmode' + str(i)][n] = interpolate(u0,V)
        obj['podmode' + str(i)][n].assign(u)
        podmode_file.close()
n_list = [2]

CU = []
for n_t in [0, 2, 4, 6]:
    print(n_t)
    CU.clear()
    for i in n_list:
    #for i in range(0,NU):
        CU_filename = 'pod_result_pre'+str(Index_power)+'/CU_'+str(i)+'/CU'+str(n_t +1)+".txt"
        Coeff = loadtxt(CU_filename)
        CU.append(Coeff)
    E_diff = 0
    i = Train_steps -1
    if i == Train_steps -1:
    #for i in range (0,Train_steps):
        a_init = Constant(0)
        a = interpolate(a_init,V)
        b = interpolate(a_init,V)
        for j in range (0, n_t +1):
            if n_t ==0:
                for n_c in range(0, len(n_list)):
                    kk = n_list[n_c]
                    a.vector().axpy(CU[n_c][i], obj['podmode' + str(kk)][j].vector())
            else:
                for n_c in range(0, len(n_list)):
                    kk = n_list[n_c]
                    a.vector().axpy(CU[n_c][i][j], obj['podmode' + str(kk)][j].vector())
        b.vector().axpy(1.9, a.vector())
        vtkfile_name = 'PVD_2d_plot_block2/Ensemble_POD_'+str(n_t +1) +'_pre'+str(Index_power)+'.pvd'
        vtkfile = File(vtkfile_name)
        vtkfile << b

