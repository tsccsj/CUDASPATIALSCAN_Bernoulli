#include <stdio.h>
#include <stdlib.h>
#define BLOCKSIZE 16
#define NumRadius 100
#define RadiusIncrement 10000

__global__ void scanKernel(float * dX, float * dY, int * dI, int nPoints, int * dWC, int * dWP, float * dLike, int nCol, int nRow, float xMin, float yMax, float cellSize, int nCase)
{
	int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
	int idInThread = threadIdx.y * blockDim.x + threadIdx.x;
	
	float cellX = xMin + cellSize * x;
	float cellY = yMax - cellSize * y;

	float dist;

	__shared__ float sX[BLOCKSIZE * BLOCKSIZE];
	__shared__ float sY[BLOCKSIZE * BLOCKSIZE];
	__shared__ int sI[BLOCKSIZE * BLOCKSIZE];

	int wCCell[NumRadius];
	int wPCell[NumRadius];

	for(int i = 0; i < NumRadius; i++)
	{
		wCCell[i] = 0;
		wPCell[i] = 0;
	}
	
	int pointProcessed;
	int pointToProcess = BLOCKSIZE * BLOCKSIZE;

	for(pointProcessed = 0; pointProcessed < nPoints; pointProcessed += BLOCKSIZE * BLOCKSIZE)
	{
		if(pointProcessed + pointToProcess > nPoints)
		{
			pointToProcess = nPoints - pointProcessed;
		}

		if(idInThread < pointToProcess)
		{
			sX[idInThread] = dX[pointProcessed + idInThread];
			sY[idInThread] = dY[pointProcessed + idInThread];
			sI[idInThread] = dI[pointProcessed + idInThread];
		}

		__syncthreads();
		for(int i = 0; i < pointToProcess; i++)
		{
			dist = sqrt((cellX - sX[i]) * (cellX - sX[i]) + (cellY - sY[i]) * (cellY - sY[i]));
			if(sI[i] > 0)
			{
				for(int j = dist / RadiusIncrement; j < NumRadius; j++)
				{
					wCCell[j] ++;
					wPCell[j] ++;
				}
			}
			else
			{
				for(int j = dist / RadiusIncrement; j < NumRadius; j++)
				{
					wPCell[j] ++;
				}
			}
		}
		__syncthreads();
	}

	int cellID = (y * nCol + x) * NumRadius;
	if(x < nCol && y < nRow)
	{
		for(int i = 0; i < NumRadius; i++)
		{
			//If it is a cold spot, ignore it (it doesn't matter to do this)
			if(wCCell[i] * nPoints < wPCell[i] * nCase)
				wCCell[i] = 0;

			dWC[cellID + i] = wCCell[i];
			dWP[cellID + i] = wPCell[i];

			dLike[cellID + i] = wCCell[i] * log((float)wCCell[i]/wPCell[i]) + (wPCell[i]-wCCell[i]) * log((float)(wPCell[i]-wCCell[i])/wPCell[i]) + (nCase-wCCell[i]) * log((float)(nCase-wCCell[i])/(nPoints-wPCell[i])) + (nPoints-wPCell[i]-nCase+wCCell[i]) * log((float)(nPoints-wPCell[i]-nCase+wCCell[i])/(nPoints-wPCell[i])); 
/*		
			if(wCCell[i] > wPCell[i])
			{
				dLike[cellID + i] = 100.0;
			}
*/
		}
	}

}

__global__ void scanKernelMC(float * dX, float * dY, int * dI, int nPoints, int * dWC, int * dAbove, int nCol, int nRow, float xMin, float yMax, int cellSize)
{
	int x = blockIdx.x * blockDim.x + threadIdx.x;
	int y = blockIdx.y * blockDim.y + threadIdx.y;
	int idInThread = threadIdx.y * blockDim.x + threadIdx.x;
	
	float cellX = xMin + cellSize * (x + 0.5);
	float cellY = yMax - cellSize * (y + 0.5);

	float dist;

	__shared__ float sX[BLOCKSIZE * BLOCKSIZE];
	__shared__ float sY[BLOCKSIZE * BLOCKSIZE];
	__shared__ int sI[BLOCKSIZE * BLOCKSIZE];

	int wCSim[NumRadius];

	for(int i = 0; i < NumRadius; i++)
	{
		wCSim[i] = 0;
	}
	
	int pointProcessed;
	int pointToProcess = BLOCKSIZE * BLOCKSIZE;

	for(pointProcessed = 0; pointProcessed < nPoints; pointProcessed += BLOCKSIZE * BLOCKSIZE)
	{
		if(pointProcessed + pointToProcess > nPoints)
		{
			pointToProcess = nPoints - pointProcessed;
		}

		if(idInThread < pointToProcess)
		{
			sX[idInThread] = dX[pointProcessed + idInThread];
			sY[idInThread] = dY[pointProcessed + idInThread];
			sI[idInThread] = dI[pointProcessed + idInThread];
		}

		__syncthreads();
		for(int i = 0; i < pointToProcess; i++)
		{
			dist = sqrt((cellX - sX[i]) * (cellX - sX[i]) + (cellY - sY[i]) * (cellY - sY[i]));
			if(sI[i] > 0)
			{
				for(int j = dist / RadiusIncrement; j < NumRadius; j++)
				{
					wCSim[j] ++;
				}
			}
		}
		__syncthreads();
	}

	int cellID = (y * nCol + x) * NumRadius;
	if(x < nCol && y < nRow)
	{
		for(int i = 0; i < NumRadius; i++)
		{
			if(dWC[cellID + i] < wCSim[i])
			{
				dAbove[cellID + i] ++;
			}
		}
	}
}


void createSample(int * ind, int nPoints, int nCase)
{
	int chosen;
	for(int i = 0; i < nPoints; i++)
		ind[i] = 0;
	while(nCase > 0)
	{
		chosen = (double)rand()/RAND_MAX * nPoints;
		if(ind[chosen] == 0)
		{
			nCase --;
			ind[chosen] = 1;
		}
	}
}

void cacScan(float * xCor, float * yCor, int * ind, int nPoints, int * wCase, int * wPop, float * like, float * pValue, int nCol, int nRow, float xMin, float yMax, float cellSize, int numSim, int nCase)
{
	float * dX;
	float * dY;
	int * dI;
	int * dWC;
	int * dWP;
	float * dLike;


	cudaError_t err;

	err = cudaMalloc((void **) &dX, nPoints * sizeof(float));
	if(err != cudaSuccess)
	{
		printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__);
		exit(1);
	}
	err = cudaMalloc((void **) &dY, nPoints * sizeof(float));
	if(err != cudaSuccess)
	{
		printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__);
		exit(1);
	}
	err = cudaMalloc((void **) &dI, nPoints * sizeof(int));
	if(err != cudaSuccess)
	{
		printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__);
		exit(1);
	}
	err = cudaMalloc((void **) &dWC, nRow * nCol * NumRadius * sizeof(int));
	if(err != cudaSuccess)
	{
		printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__);
		exit(1);
	}
	err = cudaMalloc((void **) &dWP, nRow * nCol * NumRadius * sizeof(int));
	if(err != cudaSuccess)
	{
		printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__);
		exit(1);
	}
	err = cudaMalloc((void **) &dLike, nRow * nCol * NumRadius * sizeof(float));
	if(err != cudaSuccess)
	{
		printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__);
		exit(1);
	}

	err = cudaMemcpy(dX, xCor, sizeof(float) * nPoints, cudaMemcpyHostToDevice);
	if(err != cudaSuccess)
	{
		printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__);
		exit(1);
	}
	err = cudaMemcpy(dY, yCor, sizeof(float) * nPoints, cudaMemcpyHostToDevice);
	if(err != cudaSuccess)
	{
		printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__);
		exit(1);
	}
	err = cudaMemcpy(dI, ind, sizeof(int) * nPoints, cudaMemcpyHostToDevice);
	if(err != cudaSuccess)
	{
		printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__);
		exit(1);
	}


//Kernel Goes Here
	dim3 dimBlock (BLOCKSIZE, BLOCKSIZE, 1);
	int gridX = int(ceil((float)nCol / BLOCKSIZE));
	int gridY = int(ceil((float)nRow / BLOCKSIZE));
	dim3 dimGrid (gridX, gridY, 1);

	scanKernel<<<dimGrid, dimBlock>>>(dX, dY, dI, nPoints, dWC, dWP, dLike, nCol, nRow, xMin, yMax, cellSize, nCase);

	err = cudaDeviceSynchronize();
	if(err != cudaSuccess)
	{
		printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__);
		exit(1);
	}

	err = cudaMemcpy(wCase, dWC, sizeof(int) * nCol * nRow * NumRadius, cudaMemcpyDeviceToHost);
	if(err != cudaSuccess)
	{
		printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__);
		exit(1);
	}
	err = cudaMemcpy(wPop, dWP, sizeof(int) * nCol * nRow * NumRadius, cudaMemcpyDeviceToHost);
	if(err != cudaSuccess)
	{
		printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__);
		exit(1);
	}
	err = cudaMemcpy(like, dLike, sizeof(float) * nCol * nRow * NumRadius, cudaMemcpyDeviceToHost);
	if(err != cudaSuccess)
	{
		printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__);
		exit(1);
	}


	err = cudaFree(dWP);
	if(err != cudaSuccess)
	{
		printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__);
		exit(1);
	}
	err = cudaFree(dLike);
	if(err != cudaSuccess)
	{
		printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__);
		exit(1);
	}

/*Test
	
	for(int i = 0; i < nCol * nRow * NumRadius; i++)
	{
		if(wCase[i] > wPop[i])
		{
			printf("ID: %d\twCase= %d\twPop= %d\n", i, wCase[i], wPop[i]);
		}
	}
*/

	//For simulation and pValue
	int * numAbove;
	if(NULL == (numAbove = (int *) malloc(sizeof(int) * nCol * nRow * NumRadius)))
	{
		printf("ERROR: Out of memory in line %d!\n", __LINE__);
		exit(1);
	}
	for(int i = 0; i < nCol * nRow * NumRadius; i ++)
	{
		numAbove[i] = 0;
	}

	int * dAbove;
	err = cudaMalloc((void **) &dAbove, nCol * nRow * NumRadius * sizeof(int));
	if(err != cudaSuccess)
	{
		printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__);
		exit(1);
	}
	err = cudaMemcpy(dAbove, numAbove, sizeof(int) * nCol * nRow * NumRadius, cudaMemcpyHostToDevice);
	if(err != cudaSuccess)
	{
		printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__);
		exit(1);
	}

	time_t t;
	srand((unsigned) time(&t));

	//Loop of simulation
	for(int i = 0; i < numSim; i++)
	{
		printf("Begin simulation %d\n", i);
		createSample(ind, nPoints, nCase);
		err = cudaMemcpy(dI, ind, sizeof(int) * nPoints, cudaMemcpyHostToDevice);
		if(err != cudaSuccess)
		{
			printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__);
			exit(1);
		}

		//Kernel Goes Here

		scanKernelMC<<<dimGrid, dimBlock>>>(dX, dY, dI, nPoints, dWC, dAbove, nCol, nRow, xMin, yMax, cellSize);
		
	}

	err = cudaMemcpy(numAbove, dAbove, sizeof(int) * nCol * nRow * NumRadius, cudaMemcpyDeviceToHost);
	if(err != cudaSuccess)
	{
		printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__);
		exit(1);
	}

	err = cudaFree(dWC);
	if(err != cudaSuccess)
	{
		printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__);
		exit(1);
	}
	err = cudaFree(dX);
	if(err != cudaSuccess)
	{
		printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__);
		exit(1);
	}
	err = cudaFree(dY);
	if(err != cudaSuccess)
	{
		printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__);
		exit(1);
	}
	err = cudaFree(dI);
	if(err != cudaSuccess)
	{
		printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__);
		exit(1);
	}
	err = cudaFree(dAbove);
	if(err != cudaSuccess)
	{
		printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__);
		exit(1);
	}

	for(int i = 0; i < nCol * nRow * NumRadius; i ++)
	{
		pValue[i] = (float)(numAbove[i] + 1) / (numSim + 1);
	}

	free(numAbove);
}
