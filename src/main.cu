#include <stdio.h>
#include "io.cuh"
#include "scan.cuh"

#define NumRadius 100
#define RadiusIncrement 10000

int main(int argc, char ** argv)
{
	FILE * file;

	float * xCor;
	float * yCor;
	int * ind;

	int nPoints, nCase;
	
	float xMin = -2380000, xMax = 2280000, yMin = -1520000, yMax = 1420000;
	float cellSize = 10000;

	int nRow = ceil((yMax - yMin)/cellSize);
	int nCol = ceil((xMax - xMin)/cellSize);
	yMax = yMin + nRow * cellSize;
	xMax = xMin + nCol * cellSize;

	nRow ++;
	nCol ++;

	printf("xMax = %f\txMin = %f\tyMax = %f\tyMin = %f\n", xMax, xMin, yMax, yMin);
	printf("nRow = %d\tnCol = %d\n",nRow,nCol);

	if(NULL == (file = fopen("/home/ygao29/cudaScan/data/2014_01_01", "r")))
	{
		printf("ERROR: Cannot open your flu points file \n");
		exit(1);
	}

	nPoints = getNumPoints(file);


	if(NULL == (xCor = (float *) malloc(sizeof(float)*nPoints)))
	{
		printf("ERROR: Out of memory in line %d!\n", __LINE__);
		exit(1);
	}
	if(NULL == (yCor = (float *) malloc(sizeof(float)*nPoints)))
	{
		printf("ERROR: Out of memory in line %d!\n", __LINE__);
		exit(1);
	}
	if(NULL == (ind = (int *) malloc(sizeof(int)*nPoints)))
	{
		printf("ERROR: Out of memory in line %d!\n", __LINE__);
		exit(1);
	}


	readFile(file, xCor, yCor, ind, nPoints, nCase);

	printf("nCase = %d\tnPoints = %d\n",nCase,nPoints);
	
	fclose(file);


	int * wCase;
	int * wPop;
	float * like;
	float * pValue;

	if(NULL == (wCase = (int *) malloc(sizeof(int) * nCol * nRow * NumRadius)))
	{
		printf("ERROR: Out of memory in line %d!\n", __LINE__);
		exit(1);
	}
	if(NULL == (wPop = (int *) malloc(sizeof(int) * nCol * nRow * NumRadius)))
	{
		printf("ERROR: Out of memory in line %d!\n", __LINE__);
		exit(1);
	}
	if(NULL == (like = (float *) malloc(sizeof(float) * nCol * nRow * NumRadius)))
	{
		printf("ERROR: Out of memory in line %d!\n", __LINE__);
		exit(1);
	}
	if(NULL == (pValue = (float *) malloc(sizeof(float) * nCol * nRow * NumRadius)))
	{
		printf("ERROR: Out of memory in line %d!\n", __LINE__);
		exit(1);
	}


	cacScan(xCor, yCor, ind, nPoints, wCase, wPop, like, pValue, nCol, nRow, xMin, yMax, cellSize, 99, nCase);

	free(xCor);
	free(yCor);
	free(ind);

/*
	if(NULL == (file = fopen("/home/ygao29/cudaScan/result/2014_01_01", "w")))
	{
		printf("ERROR: Cannot open your flu points file \n");
		exit(1);
	}

	for(int i = 0; i < nCol * nRow; i++)
	{
		fprintf(file, "%d:%d:%f:%f\t%d:%d:%f:%f\t%d:%d:%f:%f\t%d:%d:%f:%f\t%d:%d:%f:%f\n", wCase[i * NumRadius], wPop[i * NumRadius], like[i * NumRadius], pValue[i * NumRadius], wCase[i * NumRadius + 1], wPop[i * NumRadius + 1], like[i * NumRadius + 1], pValue[i * NumRadius + 1], wCase[i * NumRadius + 2], wPop[i * NumRadius + 2], like[i * NumRadius + 2], pValue[i * NumRadius + 2], wCase[i * NumRadius + 3], wPop[i * NumRadius + 3], like[i * NumRadius + 3], pValue[i * NumRadius + 3], wCase[i * NumRadius + 4], wPop[i * NumRadius + 4], like[i * NumRadius + 4], pValue[i * NumRadius + 4]);
	}

	fclose(file);
*/	

	int xID, yID, radiusID, ID;
	int nClusters = 50;
	float shieldFloat;
	int shieldInt;

	for(int i = 0; i < nClusters; i++)
	{
		ID = -1;
		for(int j = 0; j < nCol * nRow * NumRadius; j++)
		{
			if(wCase[j] != 0)
			{
				if(ID == -1 || like[ID] < like[j])
					ID = j;
			}
		}

		radiusID = ID % NumRadius;
		xID = (ID / NumRadius) % nCol;
		yID = ID / NumRadius / nCol;


		printf("#####################\n");
		printf("Cluster %d\n", i);
		printf("x: %f   y: %f   radius: %f\n", xMin + cellSize * xID, yMax - cellSize * yID, (float)(radiusID + 1) * RadiusIncrement);
		printf("NumPopulation: %d   NumCase: %d\n", wPop[ID], wCase[ID]);
		printf("Expected number of cases: %f\n", (float)wPop[ID] * nCase / nPoints);
		printf("Log likelihood: %f   P-Value: %f\n\n", like[ID], pValue[ID]);

		//A cluster can't be in another's center
		shieldFloat = (radiusID + 1) * RadiusIncrement / cellSize;
 		shieldInt = (int)shieldFloat;

		shieldFloat = shieldFloat * shieldFloat;

		for(int j = xID - shieldInt; j < xID + shieldInt; j++)
		{
			if(j < 0)
				j = 0;
			else if(j >= nCol)
				break;

			for(int k = yID - shieldInt; k < yID + shieldInt; k++)
			{
				if(k < 0)
					k = 0;
				else if (k >= nRow)
					break;
				if((k - yID) * (k - yID) + (j - xID) * (j - xID) < shieldFloat)
				{
					for(int l = 0; l < NumRadius; l ++)
					{
						wCase[l + NumRadius * j + NumRadius * nCol * k] = 0;
					}
				}
			}
		}
		


	}

	free(wCase);
	free(wPop);
	free(like);
	free(pValue);

	return 0;

}
