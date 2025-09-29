from fenics import *
from dolfin import *
from numpy import loadtxt
N_mode = 12
Nu = 13
for n_t in range(0,N_mode):
    E_t = 0
    E_pod = 0
    for c_m in [2, 3, 4, 6, 7, 8,9, 10, 11]:
    #for c_m in range(0,Nu):
        e_file_name = 'eigenvalue_r'+str(c_m) + '.txt'
        Eig = loadtxt(e_file_name)
        for i in range(0, len(Eig)):
            if Eig[i] > 0:
                E_t += Eig[i]
                if i > n_t:
                    E_pod += Eig[i]
    Err_file = open('theoretical_error_CPU.txt','a')
    Err_file.write('%.16g\n' % (100*sqrt(E_pod/E_t)))
    Err_file.close()
