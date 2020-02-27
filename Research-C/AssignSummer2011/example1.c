// Just trying some examples:

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

main(argc,argv)
int argc; char *argv[];
{
	int 	i,j,			a,b;
	double 			c;
	double vf;
	int entry ;	
	if(argc!=2){ 
//	printf(stderr,"\n Error %s entry", argv[0]); exit(1);}
		printf("\n NO %s entry", argv[0]);exit(1);
	}

	vf = fabs( atof(argv[1]) );
	
	printf( "\n %lf ",vf);
//	printf(" Enter the Time   :  ");
/*	if(	scanf("%d",&entry) != 1 ){
		printf(" Illegal", argv[0]);
		exit(-1);
	}
*/	entry = 5;
	for ( i=0; i<=entry; i++){
		for( j= 0; j<= entry; j++){
			printf("\n");
			printf("a");
		
		}
	}

//	printf("\n The time is : %.2d:%.2d: %3.2lf \a",a,b,c);
}
