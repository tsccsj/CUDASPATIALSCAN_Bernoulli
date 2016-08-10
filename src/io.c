#include <stdio.h>
#include <stdlib.h>

int getNumPoints(FILE * file)
{
	float x, y;
	int i;
	int count = 0;

	rewind(file);

	while(fscanf(file, "%f,%f,%d\n", &x, &y, &i) != EOF)
	{
		count ++;
	}

	return count;	
}

void readFile(FILE * file, float * xCor, float * yCor, int * ind, int nPoints, int &nCase)
{
	rewind(file);
	int indicator;

	int n = 0;
	nCase = 0;
	while(fscanf(file, "%f,%f,%d\n", (xCor + n), (yCor + n), &indicator) != EOF)
	{
		if(indicator == 1)
			nCase ++;
		ind[n] = indicator;
		n ++;
	}

	if(n != nPoints)
	{
		printf("ERROR: Incorrest input file name.\n");
		exit(1);
	}

	return;
}
