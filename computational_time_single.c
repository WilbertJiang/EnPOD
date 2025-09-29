/* import fenics C matrix and fenics G matrix and podmode. Based on podmode36*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

 
 int main(int argc,char **argv)
 { 


	double sampling_interval =4.347826086956521e-6;

        int i=0,j=0;
	double lx = 0, ux = 0.014,ly= 0,uy = 0.012,lz = 0,uz =0.0002417896,thick_actl = 0.0000557976;
	int NX = 150,NY = 150, NZ = 17,NU =50,NPOD =20,NT = 4800;
   	int NC = 20;
   	double E1,E2,E3,E4,NumError, TheoError;
	clock_t start, end;
        clock_t start_post, end_post,start_post_dev, end_post_dev;
        double cpu_time_used;
        double cpu_time_used_post,cpu_time_used_post_dev;
        double test3 = 0, test2 = 3.2, test1 = 2.5;

	/*-------------------------------------------------------------------
	load the matrix
	-------------------------------------------------------------------*/

for (int num_pod = 1; num_pod <= NPOD; num_pod++){
	 cpu_time_used = 0;
        cpu_time_used_post = 0;
        cpu_time_used_post_dev =0;


       start_post = clock();
       for(int i = 0; i< NT ; i++){
        for(int count = 0 ; count < 382500; count++){
           test3 = 0;

           for(int count1 = 0; count1 < num_pod ; count1++){
                   test3 += test2*test1;

           }
   }
}
      end_post = clock();
   cpu_time_used_post += ((double) (end_post - start_post)) / CLOCKS_PER_SEC;

        start_post_dev = clock();
	for(int i = 0; i< NT ; i++){ 
        	for(int count = 0 ; count < 67500; count++){
               		test3 = 0;
               		for(int count1 = 0; count1 < num_pod ; count1++){
                        	 test3 += test2*test1;
              	 	}
     		}
	}
    end_post_dev = clock();
    cpu_time_used_post_dev += ((double) (end_post_dev - start_post_dev)) / CLOCKS_PER_SEC;
	FILE *fileerror = fopen("Computation_time/Post_processing_time_single.txt","a");
        fprintf(fileerror,"%d\t",num_pod);
	fprintf(fileerror,"%.16lg\t", cpu_time_used_post_dev);
        fprintf(fileerror,"%.16lg\n", cpu_time_used_post);
	fclose(fileerror);

	}
	return 0;
}
	




 
