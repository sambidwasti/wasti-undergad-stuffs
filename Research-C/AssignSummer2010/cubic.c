/* 		Cubic Equation Solver               	*/
/* 	Sambid Wasti : May 26, 2010 : Gerstacker 11    	*/
/*------------------------------------------------------*/
/*	Finding the Solution to a Cubic equation	*/
/*Asking the User to input the constants in the equation*/
/*The constants are a b c and d in the equation :	*/
/*		ax^3+bx^2+cx+d=0			*/
/*	Displays the results as x1,x2 and x3		*/
/*------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
main(argc,argv)
int argc; char *argv[];
{	
	double v3,v2,v1,v0,a,b,c,d,e,f,g,h,i,j,k,l,m,n,p,r1,r2,r3,w,x,y,z;/*the variables used in the coding.*/
	printf("\nFor the equation in the form  Ax^3+Bx^2+Cx+D=0\n Input the value of the Constants A, B, C, D separated by a space: " );    
	if (scanf("%lf %lf %lf %lf",&v3,&v2,&v1,&v0)!=4)/* Specifying that only valid numbers are inputted */
	{	printf("\n%s: Improper input........ Exiting.. try again later. ",argv[0]);/* Exiting if there are some invalid input*/
		exit(-1);
	}		
	if (v3==0)/* Now inserting the quadratic coding from before */
	{	a=v2;/* matching to the old variables.*/
		b=v1;
		c=v0;	
		if ( a==0.0)	/* Specifying some cases which will reduce the time of the processor*/
		{	if (b==0.0)	
			{	printf("\n No solution for A = 0, B = 0 and C = 0 \n\n");
			}
			else if (c==0)
			{	
				printf("\n The solution is 0 \n for A = 0,B = 0 and D = 0\n\n");
			}
			else if ((b==0.0) && (c==0.0))
			{
				printf("\n There is No Equation\n\n");
			}
			else
			{	x=-c/b;
				printf("\n this is a Linear Equation \n");
				printf("\n the solution is %lf \n for A=0, B= 0,  C=%lf and D=%lf .\n \n",x,b,c);	
			}
		}
		else if (b==0.0)
		{	if (c==0.0)
			{	printf("\n the solution is 0\n for A = 0, C =0 and D=0  \n\n");
			}
			else
			{	printf("\n It is a quadratic Equation \n");
				x=-c/a;
				if(x<0)	/* to check for the complex numbers */
				{	y=sqrt(x*(-1));
					printf("\n the solutions are %lf i  and -%lf i \n for A=0,B= %lf, C=0, D=%lf\n\n",y,y,a,c);
				}
				else
				{	y=sqrt(x);
					printf("\n the solution is %lf for\n A=0,B= %lf, C=0, D=%lf\n\n",y,a,c);	
				}
			}
		}
		else if (c==0.0)
		{  	y=-b/a;
			printf("\n the solution is 0 and %lf\n for A= 0, B=%lf , C=%lf and D=0 \n\n",y,a,b);
		}
		else
		{	printf("\n It is a quadratic Equation\n ");
			z=b*b-4.0*a*c;
			if (z<0)/*checking for the complex numbers*/
			{	w=z*(-1.0);
				x=(-b/(2.0*a));
				y=((sqrt(w))/(2.0*a));
				printf("\n the solutions are %lf+%lf i    and %lf-%lf i \n for A=0, B=%lf, C=%lf and D=%lf\n\n",x,y,x,y,a,b,c);
			}
			else
			{
				w=z;
				x=(-b/(2.0*a))+((sqrt(w))/(2.0*a));
				y=(-b/(2.0*a))-((sqrt(w))/(2.0*a));
				printf("\n the solutions are %lf and %lf \n for A=0, B=%lf, C=%lf and D=%lf\n\n", x,y,a,b,c);
			}
		}
	}
	else
	{	a=v3;	/* Specifying the values in the variables */
		b=v2;
		c=v1;
		d=v0;
		f=((3.0*c/a)-((b*b)/(a*a)))/3.0;
		g=((2.0*b*b*b/(a*a*a))-(9.0*b*c/(a*a))+(27.0*d/a))/27.0;
		h=(g*g/4)+(f*f*f/27);
		/* Now we put the conditions which define the roots are real or complex */
		if ((h==0.0) && (f==0.0) && (g==0.0))    /* In this case all roots are real and equal*/
		{	r1=(cbrt(d/a))*(-1.0);
			printf("\n The solutions are all real and equal so the solution is \n X= %lf \n\n",r1);
		}
		else if (h>0) /* Case for 1 real root and others complex*/
		{	j=-(g/2.0)+(sqrt(h));
			k=cbrt(j);
			l=-((g/2.0))-sqrt(h);
			m=cbrt(l);
			r1=(k+m)-(b/(3.0*a));/* 1st Real root*/
			/* for the two complex root we have similar real parts and conjugate complex parts so*/
			n=-((k+m)/2.0)-(b/(3.0*a));
			p=((k-m)*sqrt(3.0))/2.0;
			printf("\n the solutions are \n X1=%lf \n X2= %lf + i%lf \n X3 =%lf - i%lf\n\n", r1,n,p,n,p);
		}
		else	/* Case for all real roots */
		{	i=sqrt((g*g/4.0)-h);
			j=cbrt(i);
			k=acos(-(g/(2.0*i)));
			l=j*(-1.0);
			m=cos(k/3.0);
			n=sqrt(3.0)*sin(k/3.0);
			p=(b/(3.0*a))*(-1.0);
			r1=(2.0*j*cos(k/3.0))-(b/(3.0*a));
			r2=l*(m+n)+p;
			r3=l*(m-n)+p;
			printf("\n The solutions are all real and they are \n X1= %lf \n X2= %lf \n X3= %lf\n\n",r1,r2,r3);
		}
	}
}	


