/*       Finding the Value of Pi ( The Crude Way ) 	*/
/*         Sambid Wasti: May 27,2010 (Thursday)		*/
/*------------------------------------------------------*/
/*  Specifying the two random numbers between 0 and 1   */
/*  Squaring it and addding it and comparing it with a  */
/*  circle of radius 1 and finding the value of pi via  */
/*  this method and storing the result to a file.       */
/*  We change the value of NUM  each time and get the 	*/
/*  result and store it and compare to get a better 	*/
/*  and accurate result.We also change the ISEED value  */
/*  and compare it to get a better result.		*/	
/*------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NUM 1000000000			/* number of loops we need to run so that we have that many results. 	*/
#define IMAX 10 			/* defining that the max value of a random number is 10. 		*/
#define ISEED 100000			/* Start of a random number hence called the seed and named ISEED 	*/
#define RND() (((double)random())/((double)RAND_MAX))
					/*         Here we generate the random number between 0 and 1. 		*/


main()
{
 
	int i ;				/* The variables: i for the loop 					*/
	double c,r1val, r2val, rsum, pi;/* The variables: for the two random numbers between 0 and 1.		*/
					/* and the other two for rsum and the value of pi, c is for count.	*/

	srandom (ISEED);
	c=0.0;

	for (i=0; i<NUM; i++)
	{
		r1val= RND();		/* 1st random number generated.						*/
		r2val= RND();		/* 2nd random number generated.						*/
		rsum=((r1val*r1val)+(r2val*r2val));
	
		if (rsum<=1)
		{	
			c++;		/* accepted counts gathering						*/
		}
	}

	pi=4.0*(((double)c)/((double)NUM));/* the final operation for the crude method of pi			*/
	printf("\n the value of pi = %f\n",pi);

}


  	      
