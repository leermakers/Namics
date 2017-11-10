/* File:     dot3.cu
 * Purpose:  Implement dot product on a gpu using cuda.  This version
 *           uses a binary tree reduction in which we attempt to reduce
 *           thread divergence.  It also uses shared memory to store 
 *           intermediate results.  Assumes both threads_per_block and 
 *           blocks_per_grid are powers of 2.
 *
 * Compile:  nvcc  -arch=sm_21 -o dot3 dot3.cu 
 * Run:      ./dot3 <n> <blocks> <threads_per_block>
 *              n is the vector length
 *
 * Input:    None
 * Output:   Result of dot product of a collection of random floats
 *
 */
#include <stdio.h>
#include <stdlib.h>

#define MAX_BLOCK_SZ 512

/*-------------------------------------------------------------------
 * Function:    Dev_dot  (kernel)
 * Purpose:     Implement a dot product of floating point vectors
 *              using atomic operations for the global sum
 * In args:     x, y, n
 * Out arg:     z
 *
 */
__global__ void Dev_dot(double *x, double *y, double *z, int M) {
   __shared__ double tmp[MAX_BLOCK_SZ];
   int idx = blockDim.x * blockIdx.x + threadIdx.x;
   int l_idx = threadIdx.x;
   
   if (idx < M) tmp[l_idx] = x[idx]*y[idx];
   __syncthreads();

   for (int s = blockDim.x/2; s >  0; s /= 2) {
      if (l_idx < s)
         tmp[l_idx] += tmp[l_idx + s];
      __syncthreads();
   }

   if (threadIdx.x == 0) z[blockIdx.x] = tmp[0];
   //if (l_idx == 0) z[l_idx] = tmp[0];
}  /* Dev_dot */    


/*-------------------------------------------------------------------
 * Host code 
 */
void Get_args(int argc, char* argv[], int* n_p, int* threads_per_block_p,
      int* blocks_per_grid_p);
void Setup(int n, int blocks, double** x_h_p, double** y_h_p, double** x_d_p,
      double** y_d_p, double** z_d_p);
double Serial_dot(double x[], double y[], int n);
void Free_mem(double* x_h, double* y_h, double* x_d, double* y_d,
      double* z_d);
double Dot_wrapper(double x_d[], double y_d[], double z_d[],  
      int n, int blocks, int threads);

/*-------------------------------------------------------------------
 * main
 */
int main(int argc, char* argv[]) {
   int n, threads_per_block, blocks_per_grid;
   double *x_h, *y_h, dot = 0;
   double *x_d, *y_d, *z_d;

   Get_args(argc, argv, &n, &threads_per_block, &blocks_per_grid);
   Setup(n, blocks_per_grid, &x_h, &y_h, &x_d, &y_d, &z_d);

   dot = Dot_wrapper(x_d, y_d, z_d, n, blocks_per_grid, 
         threads_per_block);

   printf("The dot product as computed by cuda is: %e\n", dot);

   dot = Serial_dot(x_h, y_h, n);
   printf("The dot product as computed by cpu is: %e\n", dot);

   Free_mem(x_h, y_h, x_d, y_d, z_d);

   return 0;
}  /* main */


/*-------------------------------------------------------------------
 * Function:  Get_args
 * Purpose:   Get and check command line args.  If there's an error
 *            quit.
 */
void Get_args(int argc, char* argv[], int* n_p, int* threads_per_block_p,
      int* blocks_per_grid_p) {

   if (argc != 4) {
      fprintf(stderr, "usage: %s <vector order> <blocks> <threads>\n", 
            argv[0]);
      exit(0);
   }
   *n_p = strtol(argv[1], NULL, 10);
   *blocks_per_grid_p = strtol(argv[2], NULL, 10);
   *threads_per_block_p = strtol(argv[3], NULL, 10);
}  /* Get_args */

/*-------------------------------------------------------------------
 * Function:  Setup
 * Purpose:   Allocate and initialize host and device memory
 */
void Setup(int n, int blocks, double** x_h_p, double** y_h_p, double** x_d_p, 
      double** y_d_p, double** z_d_p) {
   int i;
   size_t size = n*sizeof(double);

   /* Allocate input vectors in host memory */
   *x_h_p = (double*) malloc(size);
   *y_h_p = (double*) malloc(size);
   
   /* Initialize input vectors */
   srandom(1);
   for (i = 0; i < n; i++) {
      (*x_h_p)[i] = random()/((double) RAND_MAX);
      (*y_h_p)[i] = random()/((double) RAND_MAX);
   }

   /* Allocate vectors in device memory */
   cudaMalloc(x_d_p, size);
   cudaMalloc(y_d_p, size);
   cudaMalloc(z_d_p, blocks*sizeof(double));

   /* Copy vectors from host memory to device memory */
   cudaMemcpy(*x_d_p, *x_h_p, size, cudaMemcpyHostToDevice);
   cudaMemcpy(*y_d_p, *y_h_p, size, cudaMemcpyHostToDevice);
}  /* Setup */

/*-------------------------------------------------------------------
 * Function:  Dot_wrapper
 * Purpose:   CPU wrapper function for GPU dot product
 * Note:      Assumes x_d, y_d have already been 
 *            allocated and initialized on device.  Also
 *            assumes z_d has been allocated.
 */
double Dot_wrapper(double x_d[], double y_d[], double z_d[], 
      int n, int blocks, int threads) {
   int i;
   double dot = 0.0;
   double z_h[blocks];

   /* Invoke kernel */
   Dev_dot<<<blocks, threads>>>(x_d, y_d, z_d, n);
   cudaThreadSynchronize();

   cudaMemcpy(&z_h, z_d, blocks*sizeof(double), cudaMemcpyDeviceToHost);

   for (i = 0; i < blocks; i++)
      dot += z_h[i];
   return dot;
}  /* Dot_wrapper */


/*-------------------------------------------------------------------
 * Function:  Serial_dot
 * Purpose:   Compute a dot product on the cpu
 */
double Serial_dot(double x[], double y[], int n) {
   int i;
   double dot = 0;

   for (i = 0; i < n; i++)
      dot += x[i]*y[i];

   return dot;
}  /* Serial_dot */

/*-------------------------------------------------------------------
 * Function:  Free_mem
 * Purpose:   Free host and device memory
 */
void Free_mem(double* x_h, double* y_h, double* x_d, double* y_d,
      double* z_d) {

   /* Free device memory */
   cudaFree(x_d);
   cudaFree(y_d);
   cudaFree(z_d);

   /* Free host memory */
   free(x_h);
   free(y_h);

}  /* Free_mem */
