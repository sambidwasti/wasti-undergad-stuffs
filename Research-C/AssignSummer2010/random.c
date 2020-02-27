/* Random Number Generator */
/* Sambid Wasti, May 27, 2010*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Now defining the random number parameters */
#define NUM 10
#define IMAX 9
#define ISEED 110


main()
{
	int i, nval, nsum;
	double rval, rsum;
	
	srandom(ISEED);	/* Initializes the random number value generation from ISEED value that is 100 */
	
	nsum=0; rsum=0.0;
	for (i=0;i<NUM;i++)
	{	
		nval=random()%IMAX;
		printf("\n %d : \n",nval);
	}
}
