#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <algorithm> 
#include <assert.h>
#include <unistd.h>
#define DEBUG 0

typedef std::vector<int> Vector_Int;

typedef std::vector<double>Vector_Dbl;

typedef struct {

  unsigned int level;
  double xyz[6];
  Vector_Int nbr[6];
 // double *q=NULL;
	
}  cube_data;

typedef std::vector<cube_data>Cube;

/*
typedef struct 
{

  double P;
  double U;
  double V;
  double W;
	
} unknowns;

typedef std::vector<unknowns>Qs;
*/
