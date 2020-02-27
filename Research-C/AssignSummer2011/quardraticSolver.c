/* ---------------Quadratic Solver--------------------- */
/*		Sambid Wasti				*/
/*		May 16 2011 				*/
/*------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

main(argc,argv)
int argc; char *argv[];
{
	double a, b, c,
		x1, x2,
		value, svalue;	/* value for the b^2-4ac value*/

	printf( " Enter the value of the constants a, b and c ( space separated ):    ");

	if (scanf(" %lf %lf %lf" , &a, &b, &c) !=3 ){
		printf( "\n%s: Improper input ... exiting \n\n", argv[0]);
		exit(-1);
	}
	
	if ( a == 0.0 ){	 
		if ( b == 0.0 ){
			printf("\n No Solutions! ");
		}
		else {
			x1 = -c/b;
			printf("\n Solution is x = %e \n", x1);
		}
	}
	else{
		value = b*b-4*a*c;
			
		if(value <0 ){
			svalue = sqrt(-1*value);
			x1 = -b/(2*a);
			x2 = svalue/(2*a);
			printf( " \n Solutions are :\n x1 = %e + %e i \n x2 = %e - %e i\n ", x1,x2,x1,x2);	
		}
		else{
			svalue = sqrt(value);
			x1 = ( -b - svalue )/(2 *a );
			x2 = ( -b + svalue )/(2 *a );
			printf( "\n Solutions are : \n x1 = %e \n x2 = %e \n ", x1, x2);
		}
	}
}

