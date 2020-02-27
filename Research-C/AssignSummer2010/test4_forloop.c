#include <stdio.h>
#include <stdlib.h>
#include <math.h>
main()
{
	int a,b,i,j;
	int c=0;
	printf("\nPlease enter the length of the structure :\n");
	scanf("%d",&c);
	a=0;
	b=0;
	for( i =0;i<=c-1;i++)
	{
		for (j=c;j>=0;j--)
		{
			if( j>i)
			{	printf(" ");
				a++;
			}
			else
			{	printf("X");
				b++;
			}
		}
		for (j=c;j>=1;j--)
		{
			if (j>i)
			{	a++;
			}
			else
			{
			printf("X");
			b++;
			}
		}
		printf("\n");
	}
 /*	printf("%lf %lf\n",a,b);*/
}
