#!/bin/bash
for i in {0..4799}
do
	j=./solution_13blocks_pre5/file_
	k=h5
	d=./solution_13blocks_pre1_1/file_
	cp $j$i$k $d$((i+4800))$k

done
