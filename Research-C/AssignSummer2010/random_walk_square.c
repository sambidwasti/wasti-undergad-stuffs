/*----------------------------------------------*/
/*	 Random Walk on a Square Lattice 	*/
/*	 Sambid Wasti, May 27, 2010 Thursday	*/
/*----------------------------------------------*/
/* specifying 4 direction options, and using	*/
/* co-ordinate system to check the positions.	*/
/* the starting position is the origin 		*/
/*----------------------------------------------*/ 
/*   Here we are using Wasti's trick #1 where,	*/
/*   we are assuming the direction as 0,1,2,3	*/
/*   as north,south,east and west in a 2-D	*/
/*   plane.					*/		
/*----------------------------------------------*/
/* The output we need is the End-to-End Distance*/ 
/* squared and the average radius of gyration.	*/
/*----------------------------------------------*/




 			
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NUM  1000000						/* The number of loops or number of numbers generated	*/
#define IMAX  4							/* Here we only need 4 random directions.		*/
#define ISEED 999999						/* Setting up the seed as to be 1000.			*/
#define Len 2	 						/* We are setting up the length of the chain or the-	*/
								/*- no of steps taken.					*/

 

main()
{	
	srandom (ISEED);
	int i,j,d,k,l;
	double ax[Len+1],ay[Len+1],r,r1,r2,dx,dy,rssum,grsum,x,y,rsqsum,e,gyro,gyrsum,gyroadd,avggyro;
	
	rssum=0.0;
	gyroadd=0.0	;
	printf("\n");
								/* trying to take NUM  different walks of Length=LEN	*/
	for (l=0;l<NUM;l++)
	{
		x=0.0;						/* for each trial we want to start from the origin.	*/
		y=0.0;
		grsum=0.0;
		gyro=0.0;
		ax[0]=0.0;
		ay[0]=0.0;
		
		for (i=1; i<Len+1;i++)				/* for the Length or the steps of Length =Len		*/
		{	
			d=random()%IMAX;			/*Generating 4 different random directions.		*/
								/*Storing the position as a co-ordinate distances.	*/
			if (d==0)
			{	
				x++;
			}
			else if (d==1)
			{
				x--;
			}
			else if (d==2)
			{
				y++;
			}
			else
			{
				y--;
			}

				 	
			ax[i]=x;				/* Storing the Co-ordinates in an array for later use.	*/
			ay[i]=y;

		}
	
		r=((x*x)+(y*y));				/*individual distance r.				*/
		rssum+=r;				/*Summation of individual r distance squared.		*/


		/*---------------------------------------------------*/
		/* Now for the calculation of the radius of gyration.*/
		/*---------------------------------------------------*/



		for (i=0;i<Len+1;i++)				/* Defining the loop till Len				*/
		{
			
			rsqsum=0.0;				/* The counter set to 0.				*/
			
			for (j=i+1;j<Len+1;j++)			/* Similarly we did for the i values too.. here i<j	*/
			{
				dx=ax[i]-ax[j];
				dy=ay[i]-ay[j];
				r1=dx*dx+dy*dy;
				rsqsum+=r1;			/* This is the site-to-site Distance.			*/
								/* Except for dividing by N walks. We will do that later*/
			}	

			grsum+=rsqsum;				/* Taking the gyration sum where we have i<j Summation.	*/

		}			

		gyro=grsum/((((double)Len+1.0)*((double)Len+1.0)));
								/* Finding the radius of gyration for each set.		*/
		gyroadd+=gyro;					/* Adding the radius of gyrations to take the average-	*/
								/*-later.						*/
	}

	e=rssum/((double)NUM);					/* Finding the average.					*/
	printf(" For N=%i  and the length of %i we have\n Average of end to end distance square is: %lf\n",NUM,Len,e);
	
	avggyro=gyroadd/((double)NUM);				/* Finding the average radius of gyration.		*/
	printf(" And the average of the Radius of Gyration is %lf \n\n",avggyro);

}	
	
	
		








