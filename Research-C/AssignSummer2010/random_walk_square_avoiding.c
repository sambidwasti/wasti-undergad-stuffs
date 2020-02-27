/*----------------------------------------------*/
/*	 Random Walk on a Square Lattice 	*/
/*		Self Avoiding]			*/
/*	 Sambid Wasti, May 27, 2010 Thursday	*/
/*----------------------------------------------*/
/* specifying 3 direction options as it cannot 	*/
/*	take a step back., and using		*/
/* co-ordinate system to check the positions.	*/
/* the starting position is the origin 		*/
/*----------------------------------------------*/ 
/*   Here we are using Wasti's trick #1 where,	*/
/*   we are assuming the direction as 2,3,0,1	*/
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
#define ISEED 10067						/* Setting up the seed as to be 1000.			*/
#define Len 2	 						/* We are setting up the length of the chain or the-	*/
								/*- no of steps taken.					*/

 

main()
{	

	srandom (ISEED);
	int i,j,d,h,k,l,m,s,t,g;
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
		ax[0]=0.0;					/* storing location of first as origin			*/
		ay[0]=0.0;	
		s=444;						/* Any values of s except 0,1,2,3 should work.		*/
		for (i=1; i<Len+1;i++)				/* for the Length or the steps of Length =Len		*/
		{	
			g=0;
			d=random()%IMAX;			/* Generating 4 different random directions.		*/
								/* Storing the position as a co-ordinate distances.	*/
								/* Restricting the back stepping. S value assigned so	*/
								/* to compare  the old stepping. If it moved north then */
								/* next south movement should be restricted and so on...*/
			
			if((d==0) && (s==1))
			{
				g=1;				/* stored a value for g so that we can later compare.	*/
			}
			else if((d==1) && (s==0))
			{
				g=1;	
			}		
			else if((d==2) && (s==3))
			{
				g=1;
					
			}		
			else if((d==3) && (s==2))
			{
				g=1;	
				
			}		
			else if (d==0)
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
								/* Now to check if the chain comes back to its previous */
								/* Location we first check the x co-ordinate and then 	*/
								/* the y co-ordinate and see whether the co-ordinate 	*/
								/* matches or not. and if it matched then it is being 	*/
								/* repeated and that loop is cancelled by decreasing 	*/
								/* the loop value or the i value here by 1. 		*/
			

			if(g==0)				/* condion if backstepping doesnt happen.		*/
			{	
				t=0;
				for(m=0;m<i;m++)		/* Now to check if the chain comes back to its previous */
				{	
					if ((ax[m]==x) && (ay[m]==y))
					{
						t=1;		/* putting a variable so that we can check the condition*/
						break;		/* later.						*/
					}
				}
				if(t==0)			/* storing only if no conflicts are there.		*/
				{
					ax[i]=x;
					ay[i]=y;
					s=d;			/* storing the value of d in s for next comparision.	*/
				}
				else
				{
					l--;			/* if conflict is there then the walk cancelled..	*/
					break;
				}
			}	
			else					/* condition for the back-stepping.			*/
			{
				i--;
			}
		}
		if (t==0)
		{
		r=((ax[Len]*ax[Len])+(ay[Len]*ay[Len]));	/*individual distance r.				*/
		rssum+=r;					/*Summation of individual r distance squared.		*/
		
	
	

		/*---------------------------------------------------*/
		/* Now for the calculation of the radius of gyration.*/
		/*---------------------------------------------------*/



			for (i=0;i<Len+1;i++)			/* Defining the loop till Len				*/
			{
				rsqsum=0.0;			/* The counter set to 0.				*/
				for (j=i+1;j<Len+1;j++)		/* Similarly we did for the i values too.. here i<j	*/
				{
					dx=ax[i]-ax[j];
					dy=ay[i]-ay[j];
					r1=dx*dx+dy*dy;
					rsqsum+=r1;		/* This is the site-to-site Distance.			*/
								/* Except for dividing by N walks. We will do that later*/
				}	
				grsum+=rsqsum;			/* Taking the gyration sum where we have i<j Summation.	*/
			}			
			gyro=grsum/((((double)Len+1.0)*((double)Len+1.0)));
								/* Finding the radius of gyration for each set.		*/
			gyroadd+=gyro;				/* Adding the radius of gyrations to take the average-	*/
								/*-later.						*/
		}	
	}		
	e=((double)rssum)/((double)NUM);			/* Finding the average.					*/
	printf(" For N=%i  and the length of %i we have\n Average of end to end distance square is: %lf\n",NUM,Len,e);
	
	avggyro=gyroadd/((double)NUM);				/* Finding the average radius of gyration.		*/
	printf(" And the average of the Radius of Gyration is %lf \n\n",avggyro);
}
