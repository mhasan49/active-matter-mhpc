/*
File: active_matter.cu
Author: Mohammad Samsuzzaman
Institution: Department of Physics,SPPU.
Date: 2023-10-20
Description: Cuda C Code to simulate Active particles confined in a trap
Publication: https://doi.org/10.1039/C9SM01691K 
*/

#include <stdio.h>
#include <cuda_runtime.h>
#include <random>
#include <cuda.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <sys/time.h>
#include <curand.h>
#include <curand_kernel.h>
#include "./common.h"


#define MAX 100000

#define _CRT_SECURE_NO_DEPRECATE


#define num 100
#define pnum (num+1)
#define tau 100

#define a 1.0
#define bound 80.0

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

static long idum = 967369025;

float ran1(long*);
float gasdev(long*); 


__constant__ float l0 = 5.0;
__constant__ float a0_gpu=3.0;
__constant__ float a_gpu= 1.0;
__constant__ float eff = 1.0;
__constant__ float eta = 1.0;
__constant__ float fric = 15.0;
__constant__ float dt = 0.001;
__constant__ float v0 = 1.0;
__constant__ float KBT =0.00040;
__constant__ float stef =14.00;
__constant__ float hig =200.0;
__constant__ float bound_gpu =45.0;
__constant__ float alpha = 2.0;

__device__ float fpot(float r);
__device__ float f1(float r);
__device__ float x[pnum],y[pnum],vx[pnum],vy[pnum],dumx[pnum],dumy[pnum];
__device__ float direcx[pnum],direcy[pnum],vxdum[pnum],vydum[pnum],dumvx[pnum],dumvy[pnum];
__device__ float f3x=0.0f,f3y=0.0f,x2=0.0f,y2=0.0f,rr=0.0f,temp=0.0f,temp1=0.0f;
__device__ float direcabs=0.0f,r1=0.0f,r=0.0f;

/*
The function f1 takes a single argument r, which is a floating-point value representing a distance or radius.
Inside the function, r6 and r12 are defined and used for calculating powers of r. r6 represents r raised to 
the power of 6, and r12 represents r raised to the power of 12. This function calculates the force f using the 
values of eff, a_gpu, and r. 
*/

__device__ float f1(float r)
{
	float f;
	float r6;
	float r12;
	r6 = r * r * r * r * r * r;
	r12 = r6 * r6;
	f = 4 * eff * (12 * (a_gpu / r12) - 6 * (a_gpu / r6)) / r;	
	return f;
} 

/*
This function computes the trapping force as a function of the distance r from some reference point bound_gpu. 
The potential is inversely proportional to the hyperbolic cosine of the difference between r and bound_gpu.
*/

__device__ float fpot(float r)

{
	float f;
	f = -hig * stef / ((cosh(stef * (r - bound_gpu))) * (cosh(stef * (r - bound_gpu))));
	
	return f;
}



/*
The following kernel essentially performs a parallel data transfer from arrays on the GPU to other arrays on the GPU. 
The kernel is designed to be launched with a certain number of threads and blocks, which can be done using CUDA launch parameters.
*/

__global__ void position(float* x_gpu,float*y_gpu,float*vx_gpu,float*vy_gpu,int numElements )
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if ((i < numElements) && (i !=0))
	{
		
		x[i]=x_gpu[i];
		y[i]=y_gpu[i];		
		vx[i]=vx_gpu[i];
		vy[i]=vy_gpu[i];
		
		
	//	printf( "%f\t%d\n", x[i],i);
		
	}
}

__global__ void copy_position(int numElements,float *Xa, float *Ya, float *VXa,float *VYa )
{

int i = blockDim.x * blockIdx.x + threadIdx.x;
if ((i < numElements) && (i !=0))
	{	
	
		Xa[i]=x[i];
		Ya[i]=y[i];
		VXa[i]=vx[i];
		VYa[i]=vy[i];	 
}

}


/*
The following kernel updates the position and velocities of the particles.
*/

__global__ void update_position(int numElements)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < numElements)
	{	
		if(i!=0)
		{

		vx[i] = vx[i] + dumvx[i];
		vy[i] = vy[i] + dumvy[i];

		x[i] = x[i] + dumx[i];
		y[i] = y[i] + dumy[i];

		
	}
		
	} 
}

/*
The following kernel calculates the different forces acting on the particles.
*/

__global__ void vectorParent(float* gasdvx,float* gasdvy,float* rnd,int numElements )
{
	unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
	//unsigned int j = threadIdx.y + blockIdx.y * blockDim.y;

	 x2=0.0f,y2=0.0f,rr=0.0f,temp=0.0f,temp1=0.0f;
	 direcabs=0.0f,r1=0.0f,r=0.0f;

	//curandState_t state;
//	curand_init(i,0, 0, &state);
  
  
	if ((i!=0) && (i < numElements))
	{
		direcx[i] = vx[i];
		direcy[i] = vy[i];
		
			f3x = 0.0;
			f3y = 0.0;
			

	//if ((j!=0)&& (j != i)&&(j < numElements)  )
	for(int j=1;j<numElements;j++)
	{
		if(j!=i)
		{
					 x2 = x[i] - x[j];
					 y2 = y[i] - y[j];
					rr = sqrt((x2 * x2) + (y2 * y2));
					
				 	if (rr <= l0)
					{	
						direcx[i] = direcx[i] + vx[j];
						direcy[i] = direcy[i] + vy[j];

					}
					

				 	if (rr <= a0_gpu)
					{
						
						temp1 = 1 / rr;
						temp = f1(rr) * (x2 * temp1);
						f3x = f3x + temp;
						temp = f1(rr) * (y2 * temp1);
						f3y = f3y + temp;
						
					}  
 
				}
			}
		
	 	direcabs = sqrt((direcx[i] * direcx[i]) + (direcy[i] * direcy[i]));
		

//printf( "%d\n",i);
			direcx[i] = direcx[i] / direcabs;
			direcy[i] = direcy[i] / direcabs;

			r1=rnd[i];
			//r1 =  curand(&state) % MAX;
			//r1=r1/MAX;
		
			r1 = r1 - 0.50;
			r1 = r1 * eta;

			vxdum[i] = (direcx[i] * cos(r1) - direcy[i] * sin(r1));
			vydum[i] = (direcx[i] * sin(r1) + direcy[i] * cos(r1));

			r = sqrt((x[i] * x[i]) + (y[i] * y[i]));

			//r1 =  curand(&state) % MAX;
			//r1=r1/MAX;

			dumvx[i] = -fric * vx[i] * dt + v0 * vxdum[i] * dt + alpha * (fpot(r) * x[i] / r + f3x) * dt + sqrt(2.0 * fric * KBT * dt)*gasdvx[i];// * gasdev(&idum);
			dumvy[i] = -fric * vy[i] * dt + v0 * vydum[i] * dt + alpha * (fpot(r) * y[i] / r + f3y) * dt + sqrt(2.0 * fric * KBT * dt)*gasdvy[i];// * gasdev(&idum);

			dumx[i] = vx[i] * dt;
			dumy[i] = vy[i] * dt;
			 
		
			
	}

}



int
main(void)
{
	
	float L, t, dt, r, r1 = 0;
	float pi, x2, y2, theta;
	float a0;
	
	int i,j, N, nn = 0;
	N = pnum;												
	a0 = a * (pow(2.0, (1.0 / 6)));					
	L = bound - 2.0;			
	dt = 0.001;		
	t = 0.0;		
	pi = 4 * atan(1.0);	




cudaError_t err = cudaSuccess;     						
size_t nBytes = pnum * sizeof(float);

float *x_cpu = (float *)malloc(nBytes);							
float *y_cpu = (float *)malloc(nBytes);	
float *vx_cpu = (float *)malloc(nBytes);							
float *vy_cpu = (float *)malloc(nBytes);	


	i = 1;
	r1 = rand() / ((double)RAND_MAX);
	theta = 2 * pi * r1;
	r1 = rand() / ((double)RAND_MAX);
	x_cpu[i] = L * r1 * cos(theta);
	y_cpu[i] = L * r1 * sin(theta);
	r1 = rand() / ((double)RAND_MAX);
	theta = 2 * pi * r1;
	vx_cpu[i] = cos(theta);
	vy_cpu[i] = sin(theta);
	i = i + 1;
	
/*
Initialization of position and velocity
*/

while (i <= num) {
    r1 = rand() / ((double)RAND_MAX);
    theta = 2 * pi * r1;
    r1 = rand() / ((double)RAND_MAX);
    x_cpu[i] = L * r1 * cos(theta);
    y_cpu[i] = L * r1 * sin(theta);

    int collision = 0;

    for (j = 1; j <= i - 1; j++) {
        x2 = x_cpu[i] - x_cpu[j];
        y2 = y_cpu[i] - y_cpu[j];
        r = sqrt(x2 * x2 + y2 * y2);
        if (r < a0) {
            collision = 1;
            break;
        }
    }

    if (!collision) {
        r1 = rand() / ((double)RAND_MAX);
        theta = 2 * pi * r1;
        vx_cpu[i] = cos(theta);
        vy_cpu[i] = sin(theta);

        i = i + 1;
    }
}

/*  for(int i=1;i<=num;i++){
printf( "%f\t%d\n", x_cpu[i],i);
}  */

float *x_gpu = NULL;
CHECK(cudaMalloc((void **)&x_gpu, nBytes));
float *y_gpu = NULL;
CHECK(cudaMalloc((void **)&y_gpu, nBytes));
float *vx_gpu = NULL;
CHECK(cudaMalloc((void **)&vx_gpu, nBytes));
float *vy_gpu = NULL;
CHECK(cudaMalloc((void **)&vy_gpu, nBytes));

CHECK(cudaMemcpy(x_gpu, x_cpu, nBytes,  cudaMemcpyHostToDevice));
CHECK(cudaMemcpy(y_gpu, y_cpu, nBytes,  cudaMemcpyHostToDevice));
CHECK(cudaMemcpy(vx_gpu, vx_cpu, nBytes,cudaMemcpyHostToDevice));
CHECK(cudaMemcpy(vy_gpu, vy_cpu, nBytes,cudaMemcpyHostToDevice));

int threadsPerBlock = 1;
int blocksPerGrid =(pnum + threadsPerBlock - 1) / threadsPerBlock;

position<<<blocksPerGrid, threadsPerBlock>>> (x_gpu,y_gpu,vx_gpu,vy_gpu,pnum);
cudaDeviceSynchronize();
err = cudaGetLastError();

    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed vectorParent (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

float *d_x, *d_y,*d_vx, *d_vy; 
 CHECK(cudaMalloc((float**)&d_x, nBytes));
 CHECK(cudaMalloc((float**)&d_y, nBytes));
 CHECK(cudaMalloc((float**)&d_vx, nBytes));
 CHECK(cudaMalloc((float**)&d_vy, nBytes));

 float  *gpuRef;
 gpuRef  = (float *)malloc(nBytes);
 memset(gpuRef,  0, nBytes);

 CHECK(cudaMemcpy(d_x,gpuRef, nBytes, cudaMemcpyHostToDevice));
 CHECK(cudaMemcpy(d_y, gpuRef, nBytes, cudaMemcpyHostToDevice));
 CHECK(cudaMemcpy(d_vx, gpuRef, nBytes, cudaMemcpyHostToDevice));
 CHECK(cudaMemcpy(d_vy, gpuRef, nBytes, cudaMemcpyHostToDevice));
	
    int dimx = 1;
    int dimy = 1;
    dim3 block(dimx, dimy);
    dim3 grid((pnum + block.x - 1) / block.x, (pnum + block.y - 1) / block.y);
	
	//printf("grid.x %d grid.y %d grid.z %d\n", grid.x, grid.y, grid.z);
	//printf("block.x %d block.y %d block.z %d\n", block.x, block.y, block.z);
	 
 nn = 1;
 #pragma warning(suppress : 4996)
 FILE* fileOutput = fopen("/home/ubuntu/mhpc/output.txt","w");
 if (!fileOutput) {
	 perror("File opening failed");
	 return EXIT_FAILURE;
 }


 float *gasdevx_cpu,*gasdevy_cpu,*gasdevx_gpu,*gasdevy_gpu,*rand_gpu,*rand_cpu;
 gasdevx_cpu    = (float *)malloc(nBytes);
 gasdevy_cpu    = (float *)malloc(nBytes);
 rand_cpu	    = (float *)malloc(nBytes);
 
 CHECK(cudaMalloc((float**)&gasdevx_gpu, nBytes));
 CHECK(cudaMalloc((float**)&gasdevy_gpu, nBytes));
 CHECK(cudaMalloc((float**)&rand_gpu, nBytes));



 for (t = 0; t <= tau; t = t + dt)
 {
	for(int i=1;i<=num;i++)
	{
		gasdevx_cpu[i]=gasdev(&idum);
		gasdevy_cpu[i]=gasdev(&idum);
		rand_cpu[i]=rand() / ((double)RAND_MAX);
	}

	CHECK(cudaMemcpy(gasdevx_gpu, gasdevx_cpu, nBytes, cudaMemcpyHostToDevice));
	CHECK(cudaMemcpy(gasdevy_gpu, gasdevy_cpu, nBytes, cudaMemcpyHostToDevice));
	CHECK(cudaMemcpy(rand_gpu, rand_cpu, nBytes, cudaMemcpyHostToDevice));
		vectorParent<<<blocksPerGrid, threadsPerBlock>>> (gasdevx_gpu,gasdevy_gpu,rand_gpu,pnum);
	
		cudaDeviceSynchronize();

		err = cudaGetLastError();

    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed(error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

		update_position<<<blocksPerGrid, threadsPerBlock>>>(pnum);
		cudaDeviceSynchronize();
		err = cudaGetLastError();

		if (err != cudaSuccess)
		{
			fprintf(stderr, "Failed(error code %s)!\n", cudaGetErrorString(err));
			exit(EXIT_FAILURE);
		}

		if (t > (nn * 2.000))
		 {
			copy_position<<<blocksPerGrid, threadsPerBlock>>>(pnum,d_x,d_y,d_vx,d_vy);
			cudaDeviceSynchronize();

			CHECK(cudaMemcpy( x_cpu,d_x,nBytes,cudaMemcpyDeviceToHost));
			CHECK(cudaMemcpy( y_cpu,d_y,nBytes,cudaMemcpyDeviceToHost));
			CHECK(cudaMemcpy( vx_cpu,d_vx,nBytes,cudaMemcpyDeviceToHost));
			CHECK(cudaMemcpy( vy_cpu,d_vy,nBytes,cudaMemcpyDeviceToHost));
			
			
			for (int i = 1;i <= num;i++)
			{
				
				if (i == 1)
				{
					fprintf(fileOutput, "%d\n", num);
					fprintf(fileOutput, "\n");
				}
			fprintf(fileOutput, "%f\t%f\t%f\t%f\n", x_cpu[i], y_cpu[i], vx_cpu[i], vy_cpu[i]);
			
				
			}
			
			nn = nn + 1;
		}

	}
	
	fclose(fileOutput);

	/* err = cudaGetLastError();
	if (err != cudaSuccess)
	{
		printf("Failed to Launch Kernal-----");
		exit(EXIT_FAILURE);
	} */
	cudaFree(gasdevx_gpu);
	cudaFree(gasdevy_gpu);
	
	cudaFree(x_gpu);
	cudaFree(y_gpu);
	cudaFree(vx_gpu);
	cudaFree(vy_gpu);
	cudaFree(d_vx);
	cudaFree(d_vy);
	cudaFree(d_x);
	cudaFree(d_y);

	free(gasdevx_cpu);
	free(gasdevy_cpu);
	free(x_cpu);
	free(y_cpu);
	free(vx_cpu);
	free(vy_cpu);
	free(gpuRef);
	return 0;
}

/* 
The following function generates gaussianly distributed random numbers with mean 0 and variance 1
*/

float gasdev(long* idum)
{
	//float ran1(long* idum);
	static int iset = 0;
	static float gset;
	float fac, rsq, v1, v2;
	if (*idum < 0) iset = 0;
	if (iset == 0)
	{
		do {
			v1 = 2.0 * ran1(idum) - 1.0;
			v2 = 2.0 * ran1(idum) - 1.0;
			rsq = v1 * v1 + v2 * v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac = sqrt(-2.0 * log(rsq) / rsq);
		gset = v1 * fac;
		iset = 1;
		return v2 * fac;
	}
	else
	{
		iset = 0;
		return gset;
	}
}

/* 
The following function generates random number between 0 and 1.
*/

float ran1(long* idum)
{
	int j;
	long k;
	static long iy = 0;
	static long iv[32];
	float temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum = 1;
		else *idum = -(*idum);
		for (j = NTAB + 7;j >= 0;j--) {
			k = (*idum) / IQ;
			*idum = IA * (*idum - k * IQ) - IR * k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy = iv[0];
	}
	k = (*idum) / IQ;
	*idum = IA * (*idum - k * IQ) - IR * k;
	if (*idum < 0) *idum += IM;
	j = iy / NDIV;
	iy = iv[j];
	iv[j] = *idum;
	if ((temp = AM * iy) > RNMX) return RNMX;
	else return temp;
}



#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
