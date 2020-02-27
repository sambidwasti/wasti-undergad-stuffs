/* Quadratic Equation Solver*/
/* Sambid Wasti : May 25 : Gerstacker 11*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
main(argc,argv)
int argc; char *argv[];
{	
	double a,b,c,x,y,z,w,v,u;
	printf("Enter the constants A, B and C separated by a space:   ");
	if (scanf("%lf %lf %lf",&a,&b,&c)!=3)/* Specifying that only valid numbers are inputted */
	{	printf("\n%s: Improper input........ Exiting.. try again later. ",argv[0]);
			exit(-1);
	}		
	if ( a==0.0)	/* Specifying some cases which will reduce the time of the processor*/
	{	if (b==0.0)	
		{	printf("\n No solution for a = 0 and b = 0\n\n");
		}
		else if (c==0)
		{	printf("\n The solution is 0 \n for a = 0 and c = 0\n\n");
		}
		 else
		{	x=-c/b;
			printf("\n the solution is %lf \n for A=0, B=%lf and C=%lf .\n \n",x,b,c);
		}
	}
	else if (b==0.0)
	{	if (c==0.0)
		{	printf("\n the solution is 0\n for B =0 and C=0  \n\n");
		}
		else
		{	x=-c/a;
			if(x<0)	/* to check for the complex numbers */
			{	y=sqrt(x*(-1));
				printf("\n the solutions are %lf i  and -%lf i \n for A= %lf, B=0, C=%lf\n\n",y,y,a,c);
			}
			else
			{	y=sqrt(x);
				printf("\n the solution is %lf for\n A= %lf, B=0, C=%lf\n\n",y,a,c);
			}
		
		}
	}
	else if (c==0.0)
	{  	y=-b/a;
		printf("\n the solution is 0 and %lf\n for A=%lf , B=%lf and C=0 \n\n",y,a,b);
	}
	else
	{	z=b*b-4*a*c;
		if (z<0)/*checking for the complex numbers*/
		{	w=z*(-1);
			x=(-b/(2*a));
			y=((sqrt(w))/(2*a));
			printf("\n the solutions are %lf+%lf i    and %lf-%lf i \n for A=%lf, B=%lf and C=%lf\n\n",x,y,x,y,a,b,c);
		}
		else
		{
			w=z;
			x=(-b/(2*a))+((sqrt(w))/(2*a));
			y=(-b/(2*a))-((sqrt(w))/(2*a));
			printf("\n the solutions are %lf and %lf \n for A=%lf, B=%lf and C=%lf\n\n", x,y,a,b,c);
		}
	}
}



