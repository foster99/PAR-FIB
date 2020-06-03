#include "heat.h"
#include "omp.h"

/*
 * Function to copy one matrix into another
 */

void copy_mat (double *u, double *v, unsigned sizex, unsigned sizey)
{
	int blockid, howmany;
    #pragma omp parallel private(blockid,howmany)
    {
    	howmany = omp_get_num_threads();
    	blockid = omp_get_thread_num();  // Identificador del thread actual
		
		int i_start = lowerb(blockid, howmany, sizex);
      	int i_end = upperb(blockid, howmany, sizex);
		
		for (int i=max(1, i_start); i<= min(sizex-2, i_end); i++)
		    for (int j=1; j<= sizey-2; j++)
		    	v[ i*sizey+j ] = u[ i*sizey+j ];
	}
}

/*
 * Blocked Jacobi solver: one iteration step
 */
double relax_jacobi (double *u, double *utmp, unsigned sizex, unsigned sizey)
{
    double diff, sum=0.0;
    int blockid, howmany;
    #pragma omp parallel private(diff,blockid,howmany) reduction(+:sum)
    {
      howmany = omp_get_num_threads();
      blockid = omp_get_thread_num();  // Identificador del thread actual
      
      // Marcamos zona de inicio y fin para el thread iesimo (static scheduling)
      int i_start = lowerb(blockid, howmany, sizex);
      int i_end = upperb(blockid, howmany, sizex);

      for (int i=max(1, i_start); i<= min(sizex-2, i_end); i++) {
        for (int j=1; j<= sizey-2; j++) {
	     utmp[i*sizey+j]= 0.25 * ( u[ i    * sizey + (j-1) ]+  // left
	                               u[ i    * sizey + (j+1) ]+  // right
				       			   u[ (i-1)* sizey + j     ]+  // top
				       			   u[ (i+1)* sizey + j     ]); // bottom
	     diff = utmp[i*sizey+j] - u[i*sizey + j];
	     sum += diff * diff; 
	    }
      }
    }

    return sum;
}


/*
 * Blocked Gauss-Seidel solver: one iteration step
 */
double relax_gauss (double *u, unsigned sizex, unsigned sizey)
{
    double unew, diff, sum=0.0;
    int howmany = omp_get_num_threads();
    
    #pragma omp parallel for ordered(2) private(unew,diff) reduction(+:sum)
    for (int block_i=0; block_i < howmany; ++block_i){
        for(int block_j=0; block_j < howmany; ++block_j){
            
            int row_start = lowerb(block_i, howmany, sizex);
            int row_end = upperb(block_i, howmany, sizex);
            
            int col_start = lowerb(block_j, howmany, sizey);
            int col_end = upperb(block_j, howmany, sizey);
            
            #pragma omp ordered depend(sink: block_i-1,block_j) depend(sink: block_i, block_j-1)
            for (int i=max(1, row_start); i<= min(sizex-2, row_end); i++) {
                for (int j=max(1, col_start); j<= min(sizex-2, col_end); j++) {
                    unew= 0.25 * (  u[ i*sizey	+ (j-1) ]+  // left
                                    u[ i*sizey	+ (j+1) ]+  // right
                                    u[ (i-1)*sizey	+ j ]+  // top
                                    u[ (i+1)*sizey	+ j ]); // bottom
                    diff = unew - u[i*sizey+ j];
                    sum += diff * diff; 
                    u[i*sizey+j]=unew;
                }
            }
            #pragma omp ordered depend(source)
        }
    }
    return sum;
}
