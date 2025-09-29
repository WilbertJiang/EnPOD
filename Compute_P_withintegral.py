# fenics code for the thermal simulation of chip for the publication  for frequency is 3.5 GHZ
from fenics import *
from dolfin import * 
from mshr import *
from numpy import loadtxt
from petsc4py import PETSc
import numpy as np
import sys
import time
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
#define the parameter
N_STATE = 0       
Nu2 = int(sys.argv[1])# N_STATE is the controlling parameter   N_STATE = 0 indicates that we need to solve PDE and generate data; 	N_STATE = 1 indicate that we already have data just need to read it.

num_steps = 4800                                                  # num_steps is the number of time steps
Nu = 13                                                         # Nu is the number if functional unit
thick_actl =1.35e-4      
N_mode = 12
flp = loadtxt('Floorplan_AMD_multiblock.txt')
pd = loadtxt('powertrace_AMD_MLB_240_pre1_2_update.txt')                                        #pd  is power density
#compute power desity
for i in range(0,num_steps):
    for j in range(0,Nu):
        pd[i,j] = pd[i,j]/(flp[j,0]*flp[j,1]*thick_actl)
mode_inte_name = 'P_matrix_'+str(Nu2)+'.txt'
mode_inte = loadtxt(mode_inte_name)
P_matrix = np.zeros((N_mode,num_steps))
for j in range(0, N_mode):
    for i in range(0, num_steps):
        P_matrix[j][i]=pd[i,Nu2]*mode_inte[j]
header = 'P_matrix_pre1_2_update_'+str(Nu2)
P_matrix_file_name = header + '.txt'
P_matrix_file = open(P_matrix_file_name,'w')
for k1 in range(0,N_mode):
    for k2 in range(0,num_steps):
        P_matrix_file.write('%.16g\t' % (P_matrix[k1][k2]))
    P_matrix_file.write('\n')
P_matrix_file.close()

