#!/bin/bash
for i in 2 3 4 6 7 8 9 10 11
do 
	time ./ODE_solver_time $i  & 
done
#./ODE_solver_time 2 &
