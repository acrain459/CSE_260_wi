// Matrix multiply device code
#include <assert.h>
#include <math.h>
#include "utils.h"
#include "types.h"
#include <stdio.h>

#define TW 32 // = sqrt(shared memory size / number of matrices / sizeof(_DOUBLE_)) = sqrt(0xC000 / 2 / 16)
// Runs correctly for most N (but I did not test for problems)
// Not getting improved results with below results, (N = 512, 91.6 Gflops/sec) [down from ~93 Gflops]

using namespace std;

__global__ void matMul_orig(int N, _DOUBLE_ *C, _DOUBLE_ *A, _DOUBLE_ *B) 
{
	// Align Memory to 128 bits 
	/*int memsize = sizeof(_DOUBLE) * N * N;
	int alignment = 16;
	_DOUBLE_ * A_new = (_DOUBLE_ *) malloc (memsize + alignment);
	_DOUBLE_ * B_new = (_DOUBLE_ *) malloc (memsize + alignment);
	_DOUBLE_ * C_new = (_DOUBLE_ *) malloc (memsize + alignment);
	A_new += ((int) A_new % alignment) / sizeof(_DOUBLE_);
	B_new += ((int) B_new % alignment) / sizeof(_DOUBLE_);
	C_new += ((int) C_new % alignment) / sizeof(_DOUBLE_);
	memcpy (A_new, A, memsize);
	memcpy (B_new, B, memsize); //*/
	
	// Parameter Initialization
	__shared__ double As[TW][TW], Bs[TW][TW];

	int ty = threadIdx.y;
	int tx = threadIdx.x;
	int by = blockIdx.y;
	int bx = blockIdx.x;
	int I = by * TW + ty; 
	int J =  bx * TW + tx;
	
	
    //int I =  blockIdx.y*blockDim.y + threadIdx.y;
    //int J =  blockIdx.x*blockDim.x + threadIdx.x;

    if((I < N) && (J < N))
	{
		// =================================================================
        /*_DOUBLE_ _c = 0;
		
        for (unsigned int k = 0; k < N; k++) 
		{
            _DOUBLE_ a = A_new[I * N + k];
            _DOUBLE_ b = B_new[k * N + J];
            _c += a * b;
        }
		
        C_new[I * N + J] = _c;
		// ==================================================================*/
		_DOUBLE_ Cij = 0; // should be inside of loop -_-

		// TODO: round up kk loops upper bound
		// TODO: put in conditionals to avoid array out of bounds 
		for (int kk=0; kk<N/TW; kk++) // go through each block 
		{
			// read each block into shared memory
			As[ty][tx] = A[I*N + kk*TW + tx];
			Bs[ty][tx] = B[(kk*TW + ty)*N + J];
			
			__syncthreads();
			
			// dot product 
			for (int k=0; k < TW; k++)
			{
				Cij += As[ty][k] * Bs[k][tx];
			}

			__syncthreads();
		}
		
		C[I*N + J] = Cij; // <--- this does not make sense, why update every for loop without adding?  should be outside of loop 
		// ====================================================================*/
    }
	
	//memcpy (C, C_new, memsize);
}
