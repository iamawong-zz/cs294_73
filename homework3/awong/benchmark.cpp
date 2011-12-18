#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <cfloat>
#include <cmath>

#include <sys/time.h>
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
extern "C" {
#include "cblas.h"
  //#include <cblas.h>
}
#endif

/*Your function must have the following signature: */


extern void square_dgemm( int M, double *A, double *B, double *C );

/* Helper functions */

double read_timer( )
{
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = true;
    }

    gettimeofday( &end, NULL );

    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

void fill( double *p, int n )
{
    for( int i = 0; i < n; i++ )
        p[i] = 2 * drand48( ) - 1;
}

void absolute_value( double *p, int n )
{
    for( int i = 0; i < n; i++ )
        p[i] = fabs( p[i] );
}


/* The benchmarking program */

int main( int argc, char **argv )
{

    /* These sizes should highlight performance dips at multiples of certain powers-of-two */
    int test_sizes[] = {
        31, 32, 96, 97, 127, 128, 129, 191, 192, 229, 255, 256, 257,
        319, 320, 321, 417, 479, 480, 511, 512, 639, 640, 767, 768, 769,
    };
 
    /*For each test size*/
    for( unsigned int isize = 0; isize < sizeof(test_sizes)/sizeof(test_sizes[0]); isize++ )
      {
	/*Create and fill 3 random matrices A,B,C*/
        int n = test_sizes[isize];

        double *A = (double*) malloc( n * n * sizeof(double) );
        double *B = (double*) malloc( n * n * sizeof(double) );
        double *C = (double*) malloc( n * n * sizeof(double) );

        fill( A, n * n );
        fill( B, n * n );
        fill( C, n * n );
        
        /*  measure Mflop/s rate; time a sufficiently long sequence of calls to eliminate noise*/
        double Mflop_s, seconds = -1.0;
        for( int n_iterations = 1; seconds < 0.1; n_iterations *= 2 ) 
        {
            /* warm-up */
            square_dgemm( n, A, B, C );
            
            /*  measure time */
            seconds = read_timer( );
            for( int i = 0; i < n_iterations; i++ )
                square_dgemm( n, A, B, C );
            seconds = read_timer( ) - seconds;
           
            /*  compute Mflop/s rate */
            Mflop_s = 2e-6 * n_iterations * n * n * n / seconds;
        }
        printf ("Size: %d\tMflop/s: %g\n", n, Mflop_s);
        
        /*  Ensure that error does not exceed the theoretical error bound */
		
		/* Set initial C to 0 and do matrix multiply of A*B */
        memset( C, 0, sizeof( double ) * n * n );
        square_dgemm( n, A, B, C );
		/*Subtract A*B from C using standard dgemm (note that this should be 0 to within machine roundoff)*/
        cblas_dgemm( CblasColMajor,CblasNoTrans,CblasNoTrans, n,n,n, -1, A,n, B,n, 1, C,n );
		/*Subtract the maximum allowed roundoff from each element of C*/
        absolute_value( A, n * n );
        absolute_value( B, n * n );
        absolute_value( C, n * n );
        cblas_dgemm( CblasColMajor,CblasNoTrans,CblasNoTrans, n,n,n, -3.0*DBL_EPSILON*n, A,n, B,n, 1, C,n );
		/*After this test if any element in C is still positive something went wrong in square_dgemm*/
        for( int i = 0; i < n * n; i++ )
            if( C[i] > 0 )
            {
                printf( "FAILURE: error in matrix multiply exceeds an acceptable margin\n" );
                exit(-1);
            }

		/*Deallocate memory*/
        free( C );
        free( B );
        free( A );
    }
    
    return 0;
}
