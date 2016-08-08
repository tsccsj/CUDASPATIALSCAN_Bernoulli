#ifndef IOH
#define IOH

int getNumPoints(FILE * file);
void readFile(FILE * file, float * xCor, float * yCor, int * ind, int nPoints, int &nCase);

#endif
