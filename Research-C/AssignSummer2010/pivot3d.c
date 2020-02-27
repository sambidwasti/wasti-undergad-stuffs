// THE PIVOT ALGORITHM: 3-D.
// ------------------------
// Sambid Wasti. Gerstacker 11
// MON: JUNE 21 2010.
// ------------------------
// USING THE ROTATIONAL MATRIX
// TO ROTATE... 
// TAKING AN EQUILIBRIUM POINT
// AND WORKING ON IT.
// -----------------------
// BLOCK AVERAGING THE BLOCKS 
// OF 10.
// -------------------------







#include <math.h>
#include <stdio.h>
#include <stdlib.h>
					// checking floating random no.
#define NUM 1000000 
#define ISEED 10067
#define Len 4 				// Length of the chain.
#define RND() ((double)random()/(double)RAND_MAX) // For Random number between 0 and 1.

main()
{
	srandom (ISEED);
	int i,j,k,d,b,l,m,n,o,p,u,q,f,g,IMAX,accept,faccept;
	double avgend,avggyro,dx,dy,r,rx,uncgyro[10],blkend,blkgyro,ublkend,ublkgyro,bx,by,bz,cx,cy,cz,ex,ey,ez,rz,dz,diffz,avgblkend,avgblkgyro,avgublkend,avgublkgyro,uncend[10],uend,ugyro,ry,s,a,a1,a2,t,v,az[Len],ax[Len],ay[Len],ix[Len],iy[Len],iz[Len],diffx,end,endsum,uendend,ugyration,usumend,usumgyro,diffe,diffg,diffy,ang,gyroadd,gyro,rsum,grsum,rsqsum,r1;
	
				// for(i=0;i<Len;i+=radius); // this is for the different value of radius.
	j=0;
	gyroadd=0.0;


	for (i=0;i<Len;i++)	// goes from one radius to other from the origin.Now we have the initial co-ordinates.
	{
		ax[i]=j;
		ay[i]=0;
		az[i]=0;
		j++;
	}
		
	b=0;
	accept=0;
	for(k=0;k<NUM;k++)		//--------------------Equilibrium Point --------------------------/
	{	
		if(b==0)
		{	
			for(l=0;l<Len;l++)
			{
				ix[l]=ax[l];
				iy[l]=ay[l];		
				iz[l]=az[l];
			}	
		accept++;
		}
		IMAX=Len-2;		// removing the last pivot point.
		d=(random()%IMAX)+1;	// removing the first point.
		ang=RND();
		a=2*M_PI*ang;
		a1=M_PI*ang;
		a2=2*M_PI*ang;

		for(i=d+1;i<Len;i++)
				// now randomly selecting the pivot and rotating the remaining chain from that pivot.
		{
			diffx=ax[i]-ax[d];
			diffy=ay[i]-ay[d];
			diffz=az[i]-az[d];

			bx=(((cos(a2)*cos(a1)*cos(a))-(sin(a2)*sin(a)))*diffx);
			by=(((cos(a2)*cos(a1)*sin(a))+(sin(a2)*cos(a)))*diffy);
			bz=(-1*(cos(a2)*sin(a1))*diffz);
		
			cx=(((-1*sin(a2)*cos(a1)*cos(a))-(cos(a2)*sin(a)))*diffx);
			cy=(((-1*sin(a2)*cos(a1)*sin(a))+(cos(a2)*cos(a)))*diffy);
			cz=((sin(a2)*sin(a1))*diffz);
			
			ex=((sin(a1)*cos(a))*diffx);
			ey=((sin(a1)*sin(a))*diffy);
			ez=(cos(a1)*diffz);

			ax[i]=ax[d]+(bx+by+bz);
			ay[i]=ay[d]+(cx+cy+cz);
			az[i]=az[d]+(ex+ey+ez);
				
					// now checking the overlap so that it is still in the loop to reduce the no of loop
					// if there is an overlap.
			b=0;		// a counter.	
		
			for (m=0;m<d;m++)		// decrement so checking the backstepping first..
			{
							// the radius or the distance betweent he pivots.
				rx=ax[i]-ax[m];
				ry=ay[i]-ay[m];
				rz=az[i]-az[m];

				r=rx*rx+ry*ry+rz*rz;
			
				if(r<1)
				{
					b=1;
					break;
				}
			}
			
			//--- correcting or inserting the old co-ordinates if there is an overlap.	
			if(b==1)
			{
				for (n=0;n<Len;n++)
				{
					ax[n]=ix[n];
					ay[n]=iy[n];
					az[n]=iz[n];
				}
			
				break;
			}
			
			
		
			
	
		}
	}		
			//--------------------------Equilibrium Generated-------------------------//
				

	faccept=NUM/accept;
	printf("%i\n",accept);

				// For Block Averaging.
	blkend=0.0;
	blkgyro=0.0;
	ublkend=0.0;
	ublkgyro=0.0;


	for(f=0;f<10;f++)
	{
		b=0;
		endsum=0.0;
		gyroadd=0.0;
		u=0;
		for(k=0;k<NUM;k++)
		{	
			if(b==0)
			{	
				for(l=0;l<Len;l++)
				{
					ix[l]=ax[l];
					iy[l]=ay[l];
					iz[l]=az[l];
				}	
			}
			IMAX=Len-2;		// removing the last pivot point.
			d=(random()%IMAX)+1;
			ang=RND();
		
			a=2*M_PI*ang;
			a1=M_PI*ang;
			a2=2*M_PI*ang;
			
			for(i=d+1;i<Len;i++)
				// now randomly selecting the pivot and rotating the remaining chain from that pivot.
			{
				diffx=ax[i]-ax[d];
				diffy=ay[i]-ay[d];
				diffz=az[i]-az[d];

				bx=(((cos(a2)*cos(a1)*cos(a))-(sin(a2)*sin(a)))*diffx);
				by=(((cos(a2)*cos(a1)*sin(a))+(sin(a2)*cos(a)))*diffy);
				bz=(-1*(cos(a2)*sin(a1))*diffz);
			
				cx=(((-1*sin(a2)*cos(a1)*cos(a))-(cos(a2)*sin(a)))*diffx);
				cy=(((-1*sin(a2)*cos(a1)*sin(a))+(cos(a2)*cos(a)))*diffy);
				cz=((sin(a2)*sin(a1))*diffz);
				
				ex=((sin(a1)*cos(a))*diffx);
				ey=((sin(a1)*sin(a))*diffy);
				ez=(cos(a1)*diffz);

				ax[i]=ax[d]+(bx+by+bz);
				ay[i]=ay[d]+(cx+cy+cz);
				az[i]=az[d]+(ex+ey+ez);
				
			
					// now checking the overlap so that it is still in the loop to reduce the no of loop
					// if there is an overlap.
				b=0;	// a counter.	
				
				for (m=0;m<d;m++)		// decrement so checking the backstepping first..
				{
							// the radius or the distance betweent he pivots.
					rx=ax[i]-ax[m];
					ry=ay[i]-ay[m];
					rz=az[i]-az[m];
					r=rx*rx+ry*ry+rz*rz;
			
					if(r<1)
					{
						b=1;
						break;
					}
				}
				
			//--- correcting or inserting the old co-ordinates if there is an overlap.	
				if(b==1)
				{
					for (n=0;n<Len;n++)
					{
						ax[n]=ix[n];
						ay[n]=iy[n];
						az[n]=iz[n];
					}
					break;
				}
				
				
			}

		// -------- the end to end distance------- //
			{
				u++;
				end=(ax[Len-1]*ax[Len-1])+(ay[Len-1]*ay[Len-1])+(az[Len-1]*az[Len-1]);
				endsum+=end;
				uend=endsum/((double)u);
				

		
		// -------- the radius of gyration------- //
				grsum=0.0;
				for (j=0;j<Len;j++)
				{	
					rsqsum=0.0;
		
					for(o=j+1;o<Len;o++)
					{
						dx=ax[j]-ax[o];
						dy=ay[j]-ay[o];
						dz=az[j]-az[o];
						r1=dx*dx+dy*dy+dz*dz;
						rsqsum+=r1;
					}
					grsum+=rsqsum;
				}
				gyro=grsum/((((double)Len)*((double)Len)));
				gyroadd+=gyro;
				ugyro=gyroadd/((double)u);
				
			}	
		}

		avgend=endsum/((double)u);
		avggyro=gyroadd/((double)u);	
	
		uendend=sqrt(usumend/((double)NUM));
		ugyration=sqrt(usumgyro/((double)NUM));	
		blkend+=avgend;
		blkgyro+=avggyro;
		uncend[f]=blkend/((double)f+1.0);
		uncgyro[f]=blkgyro/((double)f+1.0);
	}
	avgblkend=blkend/10.0;
	avgblkgyro=blkgyro/10.0;

	for(l=0;l<10;l++)
	{
		diffe=uncend[l]-avgblkend;
		diffg=uncgyro[l]-avgblkgyro;
		usumend+=diffe*diffe;
		usumgyro+=diffg*diffg;
	}

	uendend=sqrt(usumend/((double)10.0));
	ugyration=sqrt(usumgyro/((double)10.0));	

	printf("\n For NUM: %i, ISEED: %i    and No. of Beads:%i  ",NUM,ISEED,Len); 
	printf("\n The end-to-end distance is: %lf +/- %lf\n The Radius of Gyration is : %lf +/- %lf\n",avgblkend,uendend,avgblkgyro,ugyration);
}
