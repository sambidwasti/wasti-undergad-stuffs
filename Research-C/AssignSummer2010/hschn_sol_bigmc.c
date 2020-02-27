/*--------------------------------------------------------------*/
/* File:        hschn_sol_bigmc.c				*/
/* Usage:	hschn_sol_bigmc vf D				*/
/*		(vf=solvent volume fraction, D=solvent diameter)*/
/*                                                              */
/* Monte Carlo Simulation of a Tangent Hard Sphere Chain in a	*/
/* Hard Sphere Solvent.						*/
/*								*/
/* This version is meant for large and/or high density systems. */
/* It uses a neighbor list which is updated every NUPDATE MC    */
/* cycles.  For this to be efficient, NUPDATE should be set     */
/* in range 5-20 and the resulting distance rlist must be       */
/* less than half the box size.                                 */
/*								*/
/* System consists of a single chain composed of NCHN hard 	*/
/* spheres in a solvent of NSOL hard spheres in a rectangular	*/
/* box of volume V=(Lx)x(Ly)x(Lz) with periodic boundary.	*/
/* conditions. The solvent diameter D=unit unit of length.	*/
/* vf=solvent volume fraction; rho=NSOL/V=solvent number density*/
/*								*/
/* Starting configuration: kinked chain in center of box with	*/
/* solvent particles occupying an expanded fcc lattice.		*/
/* System is equilibrated for NEQL MC cycles (1 MC cycle = N	*/
/* attempted moves) followed by NRUN cycles of data production.	*/
/* The max MC displacement is adjusted during equilibration to 	*/
/* give an acceptance fraction of about 50%.			*/
/* Data is collected in NBLK blocks for final statistics.	*/
/*								*/
/* If GFLG is set, solvent properties are computed:		*/
/* Constructs the radial distribution function g(r) and computes*/
/* the pressure (P/kTrho) by extrapolating g(r) to contact.	*/
/* g(r) is calculated out to a distance of 1+NG*DR. For accurate*/
/* pressure values a small g(r) bin width (DR) should be used.	*/
/*							7/03	*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define BLK_DMP	TRUE 	/*--write out Pij(r) functions for each block--*/
#define DIAG	TRUE	/*-- Flag for diagnostic output --*/

#define NCHN	10	/*-- Chain Length --*/
#define NSOL	400	/*-- Number of solvent particles --*/

#define NEQL    10000   /*-- # of equilibration MC cycles --*/
#define NRUN	100000	/*-- # of production MC cycles --*/
#define NBLK	10	/*-- # of blocks for statistics --*/

#define N	(NCHN+NSOL)	/*-- Total number of spheres in system --*/
#define NACCUM	5	/*-- # of MC cycles between data accumulation --*/
#define NADJUST	20	/*-- # of attempts to adjust delta during equil --*/

#define NUPDATE 10      /*-- Update date time for neighbor list --*/

/*#define ISEED	12231*/	/*-- seed for random num generator --*/

#define NBIN	10	/*--number of bins per sigma for w1N(r)--*/
#define NW	( (NCHN-2)*NBIN )  /*-total number of bins for w1N(r)-*/

#define	IVAL1	(N/4)	/*--i-value for wij(r) calc.--*/
#define	JVAL1	(3*N/4)	/*--j-value for wij(r) calc.--*/
#define NIJ1	( (JVAL1-IVAL1-1)*NBIN )

#define	IVAL2	(3*N/8)	/*--2nd i-value for wij(r) calc.--*/
#define	JVAL2	(5*N/8)	/*--2nd j-value for wij(r) calc.--*/
#define NIJ2	( (JVAL2-IVAL2-1)*NBIN )


#define GFLG	FALSE	/*-- If TRUE, compute solvent properties --*/
#define NG      150	/*-- Number of divisions in g(r) --*/
#define DR	0.020   /*-- g(r) bin width --*/

#define TRUE    1
#define FALSE   0

/*-- uniform random # in range -1.0->1.0 --*/
#define NORM    2147483647.0    /*--(2^31)-1 = largest integer--*/
#define RND()   ( (1.0-2.0*((double)random()/NORM)) )

int ncnt[NBLK];		/*--number of data points taken in each block--*/
int lst_ptr[N+1], nlist[N*N]; /*--neighbor list and pointer arrays--*/
double r2_sum[NBLK], Rg_sum[NBLK]; /*--accumulators for chain dimensions--*/
double n_sum[NBLK][NG];	/*--g(r) accumulator for each block--*/
double w1N_sum[NBLK][NW], wij1_sum[NBLK][NIJ1], wij2_sum[NBLK][NIJ2];/*--w(r) accumulator for each block--*/
double x[N], y[N], z[N]; /*--particle coordinates--*/
double Lx, Ly, Lz, rlist2, rgmax2, n_unc[NG];	/*--other global variables--*/
double X2, SIGMA, SIGMA2, SIGMA3;
double pi=M_PI;

void main(argc,argv)
int argc; char *argv[];
{
        int i, j, k, isite, nx, ny, nz, nmax, np, nblk, WIJ_FLG;
        int chn_trans(), axlrot(), chk_bnd(), overlap(), sol_ovrlp();
	double nacpt, nch1_acpt, nch2_acpt, g0, g1, gcont;
        double vp, vcp, vp_tot, rho, rho_tot, DSOL, X, V, L, delta, faccpt;
	double dth1, dth2, alpha, beta, gamma, fch1_accpt, fch2_accpt;
        double dr, davg, dnn, a, xi, yi, zi, xb, yb, zb, diff;
	double xdel, ydel, zdel, xsav,ysav,zsav,xchn[NCHN], ychn[NCHN], zchn[NCHN];
	double r, r2, r2_min, r2_max, r2_avg, r2_var, Rg_avg, Rg_var;
        double rlist, PijNRMFCT;
	double w1N[NBLK][NW], w1N_avg[NW], w1N_var[NW];
	double wij1[NBLK][NIJ1], wij1_avg[NIJ1], wij1_var[NIJ1];
	double wij2[NBLK][NIJ2], wij2_avg[NIJ2], wij2_var[NIJ2];
	void list_update(), dump_config(), g_accum(), gr_setup(), gr_output();

	/*--get volume fraction and solvent diameter from command line--*/
	if (argc!=3) {fprintf(stderr,"\nUsage: %s vf DSOL\n\n", argv[0]); exit(1);}
	vp = fabs( atof(argv[1]) );
	DSOL = atof(argv[2]);	/*--Solvent Diameter/Chain Site Diameter--*/

	SIGMA=1.0/DSOL;		/*--Chain Site Diameter--*/
	vcp=pi/(3.0*sqrt(2.0));	/*--close packed volume fraction--*/
	if (vp>vcp) {fprintf(stderr,"\nvp too large\n\n"); exit(1);}
	if (vp>0.49) printf("\nWarning: vp exceeds fluid phase stability limit\n");

	rho=(6.0/pi)*vp;	/*--particle number density--*/
	V=(double)NSOL/rho;	/*--box volume--*/
	L=cbrt(V);		/*--edge length for cubic box--*/

	X=(1.0+SIGMA)/2.0;  X2=X*X;
	SIGMA2=SIGMA*SIGMA; SIGMA3=SIGMA*SIGMA2;
	rho_tot=(double)N/V;
	vp_tot=vp+(pi/6.0)*SIGMA3*(NCHN/V);

	davg=cbrt(vcp/vp);	/*--average particle separation--*/
        delta=davg-1.0;		/*--initial max MC x,y,z displacement--*/
	dth1=pi/6.0;		/*--initial max MC angular displacement--*/
	dth2=pi/6.0;

	rlist=NUPDATE*sqrt(3.0)*delta+1.0;      /*--neighbor list radius--*/
        rlist2 = rlist*rlist;

	dr=SIGMA/(double)NBIN;
	if (IVAL1>JVAL1||IVAL1<1||JVAL1>NCHN||(JVAL1-IVAL1)==NCHN-1)      WIJ_FLG=FALSE;
	else if (IVAL2>JVAL2||IVAL2<1||JVAL2>NCHN||(JVAL2-IVAL2)==NCHN-1) WIJ_FLG=FALSE;
	else WIJ_FLG=TRUE;

	/*--Output run details--*/
        printf("\nMonte Carlo Results for %d-mer Chain\n", NCHN);
        printf("     in an N=%d Hard Sphere Fluid at solvent vp=%5.4f (vptot=%5.5f)\n",
                                                                        NSOL, vp, vp_tot);
        printf("Solvent diameter/sigma= %5.3f\n", 1.0/SIGMA);
        printf("Neighbor list in use with NUPDATE=%d\n", NUPDATE);
        printf("L=%4.2f (L/sigma=%4.2f) rho=%f (rho_tot=%f)", L, L/SIGMA, rho, rho_tot);
        if (GFLG) printf(" DR=%4.3f", DR);
        printf("\nISEED=%d NEQL=%4.2e NRUN=%4.2e NACCUM=%d\n",
                                        ISEED,(double)NEQL,(double)NRUN,NACCUM);
        printf("\ninitial delta values: dth1=%4.2f dth2=%4.2f delta=%4.2f\n", dth1, dth2, delta);
        printf("initial rlist=%3.2f\n", rlist);
        fflush(stdout);


	/*-- Setup Initial Config: fcc lattice with spacing w --*/
        dnn = 1.0+0.2*(davg-1.0);	/*--choose nn spacing w: 1<w<davg--*/
	a = sqrt(2.0)*dnn;	/*--unit cell size--*/
        nmax = (int)(L/a+0.5);
	Lx=Ly=(double)nmax*a; Lz=V/(Lx*Ly);	/*--adjust box dimensions--*/
	if (Lz>Lx) Lx=Ly=Lz=L;			/*--use cubic if possible--*/
	printf("Lx/SIG=%4.2f Ly/SIG=%4.2f Lz/SIG=%4.2f\n\n",
					Lx/SIGMA,Ly/SIGMA,Lz/SIGMA); fflush(stdout);

	/*--Place slightly kinked chain in center of box--*/
	x[NCHN/2]=Lx/2.0; y[NCHN/2]=Ly/2.0; z[NCHN/2]=Lz/2.0;
	for (i=NCHN/2-1; i>2*NCHN/6;   i--)
		{ x[i]=x[i+1]-SIGMA; y[i]=y[i+1]; z[i]=z[i+1]; }
	for (i=2*NCHN/6; i>NCHN/6;   i--)
		{ x[i]=x[i+1]; y[i]=y[i+1]-SIGMA; z[i]=z[i+1]; }
	for (i=NCHN/6; i>=0;   i--)
		{ x[i]=x[i+1]; y[i]=y[i+1]; z[i]=z[i+1]-SIGMA; }
	for (i=NCHN/2+1; i<4*NCHN/6; i++)
		{ x[i]=x[i-1]+SIGMA; y[i]=y[i-1]; z[i]=z[i-1]; }
	for (i=4*NCHN/6; i<5*NCHN/6; i++)
		{ x[i]=x[i-1]; y[i]=y[i-1]+SIGMA; z[i]=z[i-1]; }
	for (i=5*NCHN/6; i<NCHN; i++)
		{ x[i]=x[i-1]; y[i]=y[i-1]; z[i]=z[i-1]+SIGMA; }

	/*--Now place solvent particles on fcc lattice (fcc=sc+3 point basis)--*/
        for (i=NCHN, nz=0; i<N; nz++) {
	    for (ny=0; (ny<nmax && i<N); ny++) {
		for (nx=0; (nx<nmax && i<N); nx++) {
			/*--cubic lattice reference point--*/
			xb=nx*a; yb=ny*a; zb=nz*a;

			/*--try to add basis points, one at a time--*/
			x[i]=xb+0.5*a; y[i]=yb+0.5*a; z[i]=zb;
			if (overlap(i,0,NCHN)==FALSE) i++;
			if (i==N) break;

			x[i]=xb+0.5*a; y[i]=yb; z[i]=zb+0.5*a;
			if (overlap(i,0,NCHN)==FALSE) i++;
			if (i==N) break;

			x[i]=xb; y[i]=yb+0.5*a; z[i]=zb+0.5*a;
			if (overlap(i,0,NCHN)==FALSE) i++;
			if (i==N) break;

			/*--finally see if ref. point is good--*/
			x[i]=xb; y[i]=yb; z[i]=zb;
			if (overlap(i,0,NCHN)==FALSE) i++;
		}
	    }
	}
        list_update();          /*--Initial setup of neighbor list--*/
        for (i=0; i<N; i++) {	/*-- Check initial config --*/
                if (overlap(i,0,NCHN) == TRUE || sol_ovrlp(i)==TRUE) {
                        printf("Overlapping initial configuration!\n");
			dump_config(0);
                        exit(1);
                }
	}
	if (DIAG) dump_config(0);
	if (chk_bnd()==FALSE) { printf("Improper bond length\n"); fflush(stdout); }

        srandom(ISEED);			/*--initialize random num generator--*/
	if (GFLG) gr_setup(rho);	/*--Initializations for g(r) calc--*/

	/*-- SYSTEM EQULIBRATION --*/
        faccpt=fch1_accpt=fch2_accpt=0.0;
        nacpt=nch1_acpt=nch2_acpt=0.0;
        for (np=0; np<NEQL; np++) {

                if ( (np+1)%NUPDATE==0 ) list_update();

		/*-- 1 MC cycle = NCHN chn_trans() plus axlrot() attempts ... --*/
                for (i=0; i<NCHN; i++) {

			/*--Segment translation attempt--*/
                        for (j=0; j<NCHN; j++)  /*--save config--*/
                                { xchn[j]=x[j]; ychn[j]=y[j]; zchn[j]=z[j]; }
                        isite = random()%NCHN;  /*--site to be rotated--*/
                        alpha = dth1*RND();   /*--rotation angles--*/
                        beta =  (dth1/2.0)*RND();
                        gamma = dth1*RND();
                        if (chn_trans(isite,alpha,beta,gamma)==FALSE) {
                                for (j=0; j<NCHN; j++)  /*--restore config--*/
                                        { x[j]=xchn[j]; y[j]=ychn[j]; z[j]=zchn[j]; }
                        } else nch1_acpt++;       /*--or increment counter--*/

                	/*--Single site move attempt--*/
                        isite = random()%NCHN;  /*--site to be rotated--*/
			xsav=x[isite]; ysav=y[isite]; zsav=z[isite];
			if (isite==0||isite==NCHN-1) {/*--end bead rotation--*/
                        	alpha = dth2*RND();   /*--rotation angles--*/
                        	beta =  (dth2/2.0)*RND();
                        	gamma = dth2*RND();
                        	if (chn_trans(isite,alpha,beta,gamma)==FALSE) {
					x[isite]=xsav; y[isite]=ysav; z[isite]=zsav;
				} else {nch2_acpt++;}	
			} else {
				alpha = dth2*RND();
				if (axlrot(isite,alpha)==FALSE) {
					x[isite]=xsav; y[isite]=ysav; z[isite]=zsav;
				} else {nch2_acpt++;}	
			}
		}

		/*-- ... and NSOL attempted solvent moves--*/
                for (i=NCHN; i<N; i++) {
                        xi=x[i]; yi=y[i]; zi=z[i];/*--save original position--*/

			/*--Trial Move Attempt--*/
                        x[i] += delta*RND();
                        y[i] += delta*RND();
                        z[i] += delta*RND();

			/*--Check for Overlap--*/
                        if (overlap(i,0,NCHN)==TRUE
					||sol_ovrlp(i)==TRUE) { /*--REJECT: reset position--*/
                                x[i]=xi; y[i]=yi; z[i]=zi;
                        } else {	/*--ACCEPT: adjust positions for PBC's--*/
                        	if      (x[i]<0.0) x[i]+=Lx;
                        	else if (x[i]>Lx ) x[i]-=Lx;
                        	if      (y[i]<0.0) y[i]+=Ly;
                        	else if (y[i]>Ly ) y[i]-=Ly;
                        	if      (z[i]<0.0) z[i]+=Lz;
                        	else if (z[i]>Lz ) z[i]-=Lz;
				nacpt++;	/*--increment counter--*/
			}
                }
		if ((np+1)%(NEQL/NADJUST)==0) {	/*--adjust delta--*/
			faccpt = nacpt/((double)(NSOL)*(double)(NEQL/NADJUST));
			if (faccpt<0.5) delta/=1.1;
			else delta*=1.1;
                        rlist = NUPDATE*sqrt(3.0)*delta+1.0;
                        rlist2=rlist*rlist;
			nacpt=0.0;
                        printf("Equil Update: delta=%5.4f rlist=%5.4f\n",delta,rlist);
                        fflush(stdout);
		}
               	if ((np+1)%(NEQL/(NADJUST/2))==0) {     /*--adjust dth1, dth2--*/
                        fch1_accpt = nch1_acpt/(double)(NCHN*NEQL/(NADJUST/2));
                        if (fch1_accpt<0.5) dth1/=1.2;
                        else dth1*=1.2;
                        if (dth1>pi) dth1=pi;
                        fch2_accpt = nch2_acpt/(double)(NCHN*NEQL/(NADJUST/2));
                        if (fch2_accpt<0.5) dth2/=1.2;
                        else dth2*=1.2;
                        if (dth2>pi) dth2=pi;
                        nch1_acpt=nch2_acpt=0.0;
                }
        }
        printf("NEQL= %3.2e\n", (double)NEQL);
        printf("Final Equil accept: ");
        printf("ch1=%3.1f%% [dth1=%4.3f], ch2=%3.1f%% [dth2=%4.3f], ", 
        				100.0*fch1_accpt, dth1, 100.0*fch2_accpt, dth2);
        printf("fld=%3.1f%% [delta=%4.3f]  RLIST=%5.3f\n\n", 100.0*faccpt, delta, rlist);
	fflush(stdout);

	if (chk_bnd()==FALSE) printf("Improper bond length\n");
	if (DIAG) dump_config(NEQL);

	/*-- DATA PRODUCTION --*/
        nacpt=nch1_acpt=nch2_acpt=0.0; r2_min=(NCHN-1)*(NCHN-1)*SIGMA2; r2_max=0.0;
        for (np=0; np<NRUN; np++) {	/*--one MC cycle = N attempted moves--*/

               if ( np%NUPDATE==0 ) list_update();

		/*--CHAIN MOVES--*/
                for (i=0; i<NCHN; i++) {

			/*--Segment translation attempt--*/
                        for (j=0; j<NCHN; j++)  /*--save config--*/
                                { xchn[j]=x[j]; ychn[j]=y[j]; zchn[j]=z[j]; }
                        isite = random()%NCHN;  /*--site to be rotated--*/
                        alpha = dth1*RND();   /*--rotation angles--*/
                        beta =  (dth1/2.0)*RND();
                        gamma = dth1*RND();
                        if (chn_trans(isite,alpha,beta,gamma)==FALSE)
                                for (j=0; j<NCHN; j++)  /*--restore config--*/
                                        { x[j]=xchn[j]; y[j]=ychn[j]; z[j]=zchn[j]; }
                        else nch1_acpt++;       /*--or increment counter--*/

                 	/*--Single site move attempts--*/
                        isite = random()%NCHN;  /*--site to be rotated--*/
			xsav=x[isite]; ysav=y[isite]; zsav=z[isite];
			if (isite==0||isite==NCHN-1) {
                        	alpha = dth2*RND();   /*--rotation angles--*/
                        	beta =  (dth2/2.0)*RND();
                        	gamma = dth2*RND();
                        	if (chn_trans(isite,alpha,beta,gamma)==FALSE) {
					x[isite]=xsav; y[isite]=ysav; z[isite]=zsav;
				} else {nch2_acpt++;}	
			} else {
				alpha = dth2*RND();
				if (axlrot(isite,alpha)==FALSE) {
					x[isite]=xsav; y[isite]=ysav; z[isite]=zsav;
				} else {nch2_acpt++;}	
			}
		}

		/*--FLUID MOVES--*/
                for (i=NCHN; i<N; i++) {/*--select disks in order--*/
                        xi=x[i]; yi=y[i]; zi=z[i]; /*--save original position--*/

			/*--Trial Move Attempt--*/
                        x[i] += delta*RND();
                        y[i] += delta*RND();
                        z[i] += delta*RND();

			/*--Check for Overlap--*/
                        if (overlap(i,0,NCHN)==TRUE
					||sol_ovrlp(i)==TRUE) { /*--REJECT: reset position--*/
                                x[i]=xi; y[i]=yi; z[i]=zi;
                        } else {	/*--ACCEPT: adjust positions for pbc--*/
                        	if      (x[i]<0.0) x[i]+=Lx;
                        	else if (x[i]>Lx ) x[i]-=Lx;
                        	if      (y[i]<0.0) y[i]+=Ly;
                        	else if (y[i]>Ly ) y[i]-=Ly;
                        	if      (z[i]<0.0) z[i]+=Lz;
                        	else if (z[i]>Lz ) z[i]-=Lz;
				nacpt++;	/*--increment counter--*/
			}
                }

                if (np%NACCUM==0) {      /*--take data every NACCUM cycle-*/
                        nblk=np/(NRUN/NBLK);    /*--current block number--*/
			if (nblk>=NBLK) {printf("Error: nblk too large!\n"); fflush(stdout); }
                        ncnt[nblk]++;           /*--increment counter--*/

			if (GFLG) g_accum(nblk);	/*--accumulate data--*/

                        xdel=x[NCHN-1]-x[0]; ydel=y[NCHN-1]-y[0]; zdel=z[NCHN-1]-z[0];
                        r2 = (xdel*xdel+ydel*ydel+zdel*zdel);
			if (r2<r2_min) r2_min=r2;
			if (r2>r2_max) r2_max=r2;
                        r2_sum[nblk] += r2;     /*--end-to-end distance--*/

			/*--build histogram for P1N(r)--*/
			k=(int)((sqrt(r2)-SIGMA)/dr);
			if (k<NW) w1N_sum[nblk][k]++;

			/*--build histogram for Pij(r)--*/
			if (WIJ_FLG) {
				i=IVAL1-1; j=JVAL1-1;	/*--1st wij ... i<j contribution--*/
                        	xdel=x[j]-x[i]; ydel=y[j]-y[i]; zdel=z[j]-z[i];
                        	r2 = (xdel*xdel+ydel*ydel+zdel*zdel);
				k=(int)((sqrt(r2)-SIGMA)/dr);
				if (k<NIJ1) wij1_sum[nblk][k]++;

				i=NCHN-JVAL1; j=NCHN-IVAL1; /*--i>j contribution--*/
                        	xdel=x[j]-x[i]; ydel=y[j]-y[i]; zdel=z[j]-z[i];
                        	r2 = (xdel*xdel+ydel*ydel+zdel*zdel);
				k=(int)((sqrt(r2)-SIGMA)/dr);
				if (k<NIJ1) wij1_sum[nblk][k]++;

				i=IVAL2-1; j=JVAL2-1;	/*--2nd wij ... i<j contribution--*/
                        	xdel=x[j]-x[i]; ydel=y[j]-y[i]; zdel=z[j]-z[i];
                        	r2 = (xdel*xdel+ydel*ydel+zdel*zdel);
				k=(int)((sqrt(r2)-SIGMA)/dr);
				if (k<NIJ2) wij2_sum[nblk][k]++;

				i=NCHN-JVAL2; j=NCHN-IVAL2; /*--i>j contribution--*/
                        	xdel=x[j]-x[i]; ydel=y[j]-y[i]; zdel=z[j]-z[i];
                        	r2 = (xdel*xdel+ydel*ydel+zdel*zdel);
				k=(int)((sqrt(r2)-SIGMA)/dr);
				if (k<NIJ2) wij2_sum[nblk][k]++;
			}

                        for (i=0; i<NCHN-1; i++)
                                for (j=i+1; j<NCHN; j++) {/*-radius of gyration-*/
                                        xdel=x[j]-x[i]; ydel=y[j]-y[i]; zdel=z[j]-z[i];
                                        r2 = (xdel*xdel+ydel*ydel+zdel*zdel);
                                        Rg_sum[nblk] += r2/(double)(NCHN*NCHN);
                                }
                        /*--Output data after each block complete--*/
                        if (ncnt[nblk]==(NRUN/(NBLK*NACCUM))) {
                                printf("Block %2d complete: ", nblk+1);
                                printf("<r2>= %6.4f  <Rg2>= %6.4f     ",
                                        	r2_sum[nblk]/(SIGMA2*ncnt[nblk]),
							Rg_sum[nblk]/(SIGMA2*ncnt[nblk]));
				printf("[r2min= %6.4f  r2max= %6.4f]\n",
							r2_min/SIGMA2, r2_max/SIGMA2);
				if (GFLG) {
                		    g0 = n_sum[nblk][0]/(n_unc[0]*(double)(NSOL*ncnt[nblk]));
                		    g1 = n_sum[nblk][1]/(n_unc[1]*(double)(NSOL*ncnt[nblk]));
        			    gcont = 1.5*g0-0.5*g1; /*-linear extrap for g_contact-*/
                                    printf("solvent pressure= %6.4f\n", 1.0+4.0*vp*gcont);
				}
				if (BLK_DMP) { /*--write out Pij(r) functions for this block--*/
                			PijNRMFCT = (SIGMA/dr)/(ncnt[nblk]);
        				for (i=0; i<NW; i++) {
					    r=SIGMA+((double)i+0.5)*dr;
					    /*--construct P1N(r)--*/
                			    printf("%5.4f  %6.4f", r/SIGMA, PijNRMFCT*w1N_sum[nblk][i]);
					    if (i<NIJ1)  /*--construct Pij1(r)--*/
                				printf("   %6.4f", PijNRMFCT*wij1_sum[nblk][i]/2.0);
					    if (i<NIJ2)  /*--construct Pij213(r)--*/
                				printf("   %6.4f", PijNRMFCT*wij2_sum[nblk][i]/2.0);
					    printf("\n");
					}
					printf("\n");
				}
                                fflush(stdout);
				/*--reset r2 min, max --*/
        			r2_min=(NCHN-1)*(NCHN-1)*SIGMA2; r2_max=0.0;
                        }
		}


        }
	/*-- DATA PRODUCTION COMPLETE --*/
	
	list_update();
        for (i=0; i<N; i++)	/*--Check final configuration for overlaps--*/
                if (overlap(i,0,N)==TRUE) {
                        printf("Overlapping final configuration!\n");
                        exit(1);
                }
	if (chk_bnd()==FALSE) printf("Improper bond length\n");

        faccpt = nacpt/((double)(NSOL)*(double)(NRUN));
        fch1_accpt = nch1_acpt/((double)(NCHN)*(double)(NRUN));
        fch2_accpt = nch2_acpt/((double)(NCHN)*(double)(NRUN));
        printf("NRUN= %3.2e  (Final Prod accept: ", (double)NRUN);
        printf("ch1=%3.1f%%, ch2=%3.1f%%, fld=%3.1f%%)\n\n", 
				100.0*fch1_accpt, 100.0*fch2_accpt, 100.0*faccpt);

	if (DIAG) dump_config(NEQL+NRUN);

	/*---CONSTRUCT AND OUTPUT FINAL RESULTS---*/
	if (GFLG) gr_output(vp);

        r2_avg=Rg_avg=0.0;
        for (i=0; i<NW; i++) {
		w1N_avg[i]=0.0;
		if (WIJ_FLG && i<NIJ1) wij1_avg[i]=0.0;
		if (WIJ_FLG && i<NIJ2) wij2_avg[i]=0.0;
	}
        for (nblk=0; nblk<NBLK; nblk++) {
                r2_avg += r2_sum[nblk]/ncnt[nblk];
                Rg_avg += Rg_sum[nblk]/ncnt[nblk];
        	for (i=0; i<NW; i++) {	/*--construct P1N(r)--*/
                	w1N[nblk][i] = (SIGMA/dr)*w1N_sum[nblk][i]/(ncnt[nblk]);
			w1N_avg[i] += w1N[nblk][i];
			if (WIJ_FLG && i<NIJ1) {	/*--construct 1st Pij(r)--*/
                		wij1[nblk][i] = (SIGMA/dr)*wij1_sum[nblk][i]/(2.0*ncnt[nblk]);
				wij1_avg[i] += wij1[nblk][i];
			}
			if (WIJ_FLG && i<NIJ2) {	/*--construct 2nd Pij(r)--*/
                		wij2[nblk][i] = (SIGMA/dr)*wij2_sum[nblk][i]/(2.0*ncnt[nblk]);
				wij2_avg[i] += wij2[nblk][i];
			}
		}
        }
        r2_avg/=NBLK; Rg_avg/=NBLK;
        for (i=0; i<NW; i++) {
		w1N_avg[i] /= NBLK;
		if (WIJ_FLG && i<NIJ1) wij1_avg[i] /= NBLK;
		if (WIJ_FLG && i<NIJ2) wij2_avg[i] /= NBLK;
	}

        r2_var=Rg_var=0.0;
        for (i=0; i<NW; i++) {
		w1N_var[i]=0.0;
		if (WIJ_FLG && i<NIJ1) wij1_var[i]=0.0;
		if (WIJ_FLG && i<NIJ2) wij2_var[i]=0.0;
	}
        for (nblk=0; nblk<NBLK; nblk++) {
                diff = r2_avg-r2_sum[nblk]/ncnt[nblk];
                r2_var += diff*diff;
                diff = Rg_avg-Rg_sum[nblk]/ncnt[nblk];
                Rg_var += diff*diff;
        	for (i=0; i<NW; i++) {
			diff = w1N_avg[i]-w1N[nblk][i];
			w1N_var[i] += diff*diff;
			if (WIJ_FLG && i<NIJ1) {
				diff = wij1_avg[i]-wij1[nblk][i];
				wij1_var[i] += diff*diff;
			}
			if (WIJ_FLG && i<NIJ2) {
				diff = wij2_avg[i]-wij2[nblk][i];
				wij2_var[i] += diff*diff;
			}
		}
        }
        r2_var/=(NBLK-1); Rg_var/=(NBLK-1);
        for (i=0; i<NW; i++) {
		w1N_var[i]/=(NBLK-1);
		if (WIJ_FLG && i<NIJ1) wij1_var[i]/=(NBLK-1);
		if (WIJ_FLG && i<NIJ2) wij2_var[i]/=(NBLK-1);
	}

	r2_avg/=SIGMA2; r2_var/=(SIGMA2*SIGMA2);
	Rg_avg/=SIGMA2; Rg_var/=(SIGMA2*SIGMA2);
        printf("<r1n^2> = %6.4f +- %5.4f\n", r2_avg, sqrt(r2_var));
        printf("<Rg^2>  = %6.4f +- %5.4f\n", Rg_avg, sqrt(Rg_var));

	if (WIJ_FLG) printf("\n  r       P1N(r)               P%d%d(r)               P%d%d(r)\n",
									IVAL1,JVAL1,IVAL2,JVAL2);
	else 	     printf("\n  r       P1N(r)\n");
        for (i=0; i<NW; i++) {	/*--r-value for w(r)=center of w[i] bin--*/
		r=SIGMA+((double)i+0.5)*dr;
                printf("%5.4f  %5.3f +/- %5.3f    ", r/SIGMA, w1N_avg[i], sqrt(w1N_var[i]));
                if (WIJ_FLG&&i<NIJ1) printf("  %5.3f +/- %5.3f    ", wij1_avg[i], sqrt(wij1_var[i]));
                if (WIJ_FLG&&i<NIJ2) printf("  %5.3f +/- %5.3f    ", wij2_avg[i], sqrt(wij2_var[i]));
		printf("\n");
	}
}

int overlap(i,nstrt,nstop) /*--check for overlap of particle i--*/
int i,nstrt,nstop;	   /*--PBC's enforced even if chain outside the box--*/
{
        int j;
        double xi,yi,zi,rx,ry,rz,r2,d2;

        xi=x[i]; yi=y[i]; zi=z[i];
        for (j=nstrt; j<nstop; j++) {

                if (j==i) continue;
		/*--set d2min for current ij pair--*/
		if (i<NCHN&&j<NCHN) {
			if (abs(i-j)==1) continue; /*--skip on bonded sites--*/
			else d2=SIGMA2;			/*--chain-chain--*/
		}
		else if (i>=NCHN&&j>=NCHN) d2=1.0;	/*--solvent-solvent--*/
		else			   d2=X2;	/*--solvent-chain--*/

                rx = fabs(xi-x[j]);
                ry = fabs(yi-y[j]);
                rz = fabs(zi-z[j]);
                while (rx>Lx/2.0) rx -= Lx;
                while (ry>Ly/2.0) ry -= Ly;
                while (rz>Lz/2.0) rz -= Lz;
                r2 = rx*rx + ry*ry + rz*rz;
                if (r2<d2) return(TRUE);	/*--return if overlap--*/
        }
        return(FALSE);	/*--No overlap if you make it here--*/
}

int sol_ovrlp(i)	/*--check for sol-sol overlaps using nlist--*/
int i;
{
        int j,jl,jlstrt,jlstop;
        double xi,yi,zi,rx,ry,rz,r2;

        xi=x[i]; yi=y[i]; zi=z[i];

        jlstrt=lst_ptr[i]; jlstop=lst_ptr[i+1];
        for (jl=jlstrt; jl<jlstop; jl++) {
                j = nlist[jl];
                if (j==i) continue;
                rx = fabs(xi-x[j]);
                ry = fabs(yi-y[j]);
                rz = fabs(zi-z[j]);
                while (rx>Lx/2.0) rx -= Lx;
                while (ry>Ly/2.0) ry -= Ly;
                while (rz>Lz/2.0) rz -= Lz;
                r2 = rx*rx + ry*ry + rz*rz;
                if (r2<1.0) return(TRUE);
        }
        return(FALSE);
}

void list_update()
{
        int i,j,ptr_val;
        double xi,yi,zi,rx,ry,rz,r2;

        ptr_val=0;
        for (i=NCHN; i<N; i++) {
                lst_ptr[i]=ptr_val;
                xi=x[i]; yi=y[i]; zi=z[i];
                for (j=NCHN; j<N; j++) {
                        if (j==i) continue;
                        rx = fabs(xi-x[j]);
                        ry = fabs(yi-y[j]);
                        rz = fabs(zi-z[j]);
                        while (rx>Lx/2.0) rx -= Lx;
                        while (ry>Ly/2.0) ry -= Ly;
                        while (rz>Lz/2.0) rz -= Lz;
                        r2 = rx*rx + ry*ry + rz*rz;
                        if (r2 <= rlist2) { nlist[ptr_val]=j; ptr_val++; }
                        if (r2<1.0) {printf("overlap in list_update()!!\n"); fflush(stdout);}
                }
        }
        lst_ptr[N]=ptr_val;
}


int chk_bnd()	/*--check chain bondlengths--*/
{
	int i;
	double xdel, ydel, zdel, r2;

	for (i=1; i<NCHN; i++) {
		xdel=x[i]-x[i-1]; ydel=y[i]-y[i-1]; zdel=z[i]-z[i-1];
		r2 = xdel*xdel + ydel*ydel + zdel*zdel;
		if (fabs(r2-SIGMA2)>1.0e-8) return(FALSE);
	}
	return(TRUE);
}

int chn_trans(is,alpha,beta,gamma)   /*--translation move for site is sub-chain--*/
int is; double alpha,beta,gamma;
{
        int j, ip, overlap();
        double xi,yi,zi,xip,yip,zip,xdel,ydel,zdel;
	double sin_a, cos_a, sin_b, cos_b, sin_g, cos_g;
        double Rxx,Rxy,Rxz,Ryx,Ryy,Ryz,Rzx,Rzy,Rzz;

        xi=x[is]; yi=y[is]; zi=z[is];
        if (is<NCHN/2)  ip=is+1;
        else            ip=is-1;
        xip=x[ip]; yip=y[ip]; zip=z[ip];

        /*--Get displacement associated with site is rotation--*/
        sin_a=sin(alpha); sin_b=sin(beta); sin_g=sin(gamma);
        cos_a=cos(alpha); cos_b=cos(beta); cos_g=cos(gamma);
        Rxx = (cos_g*cos_b*cos_a-sin_g*sin_a);
        Rxy = (cos_g*cos_b*sin_a+sin_g*cos_a);
        Rxz = (-cos_g*sin_b);
        Ryx = (-sin_g*cos_b*cos_a-cos_g*sin_a);
        Ryy = (-sin_g*cos_b*sin_a+cos_g*cos_a);
        Ryz = (sin_g*sin_b);
        Rzx = (sin_b*cos_a);
        Rzy = (sin_b*sin_a);
        Rzz = (cos_b);

        xdel = xip-xi + (xi-xip)*Rxx + (yi-yip)*Rxy + (zi-zip)*Rxz;
        ydel = yip-yi + (xi-xip)*Ryx + (yi-yip)*Ryy + (zi-zip)*Ryz;
        zdel = zip-zi + (xi-xip)*Rzx + (yi-yip)*Rzy + (zi-zip)*Rzz;

        /*--Translate short end of chain, checking for intrachain overlaps--*/
        if (is<NCHN/2) { /*--move site 0->is end--*/
                for (j=is; j>=0; j--) {
                        x[j]+=xdel; y[j]+=ydel; z[j]+=zdel;
                        if (overlap(j,ip,NCHN)==TRUE) return(FALSE);
                }
        } else { /*--move site is->NCHN end--*/
                for (j=is; j<NCHN; j++) {
                        x[j]+=xdel; y[j]+=ydel; z[j]+=zdel;
                        if (overlap(j,0,is)==TRUE) return(FALSE);
                }
        }

        /*--Finally check for intermolecular overlap--*/
        if (is<NCHN/2) {  /*--check either for chain sites 0->is ...--*/
                for (j=is; j>=0; j--)
                        if (overlap(j,NCHN,N)==TRUE) return(FALSE);
        } else {                /*--... or for chain sites is->NCHN--*/
                for (j=is; j<NCHN; j++)
                        if (overlap(j,NCHN,N)==TRUE) return(FALSE);
        }

        return(TRUE); /*--No overlaps if you get to here--*/
}

int axlrot(is,alpha)     /*--axial rotation of site is--*/
int is; double alpha;
{
        int i,j,k,l;
        int overlap();
        double xp1,xp2,yp1,yp2,zp1,zp2,xi,yi,zi;
        double xdel,ydel,zdel,r,s,cos_p,sin_p,cos_t,sin_t,cos_a,sin_a;
        double a[3][3],b[3][3],c[3][3];

        cos_a=cos(alpha); sin_a=sin(alpha);

        /*--Determine (theta,phi) orientation of the rotation axis--*/
        xp1=x[is-1]; yp1=y[is-1]; zp1=z[is-1];
        xp2=x[is+1]; yp2=y[is+1]; zp2=z[is+1];
        xdel=xp2-xp1; ydel=yp2-yp1; zdel=zp2-zp1;
        r=sqrt(xdel*xdel+ydel*ydel+zdel*zdel);
        s=sqrt(xdel*xdel+ydel*ydel);

	if (s!=0.0) {
        	cos_t=xdel/s; sin_t=ydel/s; /*--0<theta<2pi (the azimuth angle)--*/
        	cos_p=zdel/r; sin_p=s/r;    /*-- 0<phi<pi (the polar angle) --*/
	} else {
		cos_t=1.0; sin_t=1.0; cos_p=1.0; sin_p=0.0;	/*--s=0.0--*/
	}

        /*--Setup the alpha rotation Matrix A--*/
        a[0][0]=cos_a;    a[0][1]=sin_a; a[0][2]=0.0;
        a[1][0]=(-sin_a); a[1][1]=cos_a; a[1][2]=0.0;
        a[2][0]=0.0;      a[2][1]=0.0;   a[2][2]=1.0;

        /*--Setup the (theta,phi) rotation Matrix B--*/
        b[0][0]=cos_p*cos_t; b[0][1]=cos_p*sin_t; b[0][2]=(-sin_p);
        b[1][0]=(-sin_t);    b[1][1]=cos_t;       b[1][2]=0.0;
        b[2][0]=sin_p*cos_t; b[2][1]=sin_p*sin_t; b[2][2]=cos_p;

        /*--Setup the full rotation Matrix C=(B^T)A(B)--*/
        for (i=0; i<3; i++)
            for (j=0; j<3; j++) {
                c[i][j]=0.0;
                for (k=0; k<2; k++)
                    for (l=0; l<2; l++)
                        c[i][j] += b[k][i]*a[k][l]*b[l][j];
                c[i][j] += b[2][i]*b[2][j];
            }

        /*--rotate site is about ip1-ip2 axis--*/
        xi=x[is]-xp1; yi=y[is]-yp1; zi=z[is]-zp1;
        x[is] = xp1+c[0][0]*xi+c[0][1]*yi+c[0][2]*zi;
        y[is] = yp1+c[1][0]*xi+c[1][1]*yi+c[1][2]*zi;
        z[is] = zp1+c[2][0]*xi+c[2][1]*yi+c[2][2]*zi;

        if (overlap(is,0,N) == TRUE) return(FALSE);
	else return(TRUE);
}

void dump_config(ncycle)	/*--Output current system configuration--*/
int ncycle;
{
	int i;

	printf("\nSystem Config after %d MC cycles:\n", ncycle);
	for (i=0; i<N; i++)
		printf("%d %f %f %f\n", i, x[i], y[i], z[i]); 
	printf("\n"); fflush(stdout);
}



void gr_setup(rho) /*-- Some initializations for g(r) construction --*/
double rho;
{
	int i;
	double nsol, ai, rgmax, rmin, rmax, Vshell;

	nsol = (double)NSOL;
        for (i=0; i<NG; i++) {	/*--normalization for each g(r) shell--*/
                ai = (double)i;
                rmin = 1.0+ai*DR;
                rmax = rmin+DR;
		Vshell = 4.0*pi*(rmax*rmax*rmax-rmin*rmin*rmin)/3.0;
                n_unc[i] = Vshell*rho*(nsol-1.0)/nsol; /*--uncorrelated fluid result--*/
        }
	/*--Set max r value for g(r) ... should not exceed half box width--*/
        rgmax=1.0+(double)NG*DR;
	if (rgmax>Lx/2.0) rgmax=Lx/2.0;
	if (rgmax>Ly/2.0) rgmax=Ly/2.0;
	if (rgmax>Lz/2.0) rgmax=Lz/2.0;
        rgmax2=rgmax*rgmax;	/*--global variable used in gaccum()--*/
}

void g_accum(nblk)	/*--g(r) bin accumulator--*/
int nblk;
{
        int i,j,rindx;
        double xi,yi,zi,rx,ry,rz,r,r2;

        for (i=NCHN; i<N; i++) {   /*--sweep through all NSOL particles ...--*/
        	xi=x[i]; yi=y[i]; zi=z[i]; /*--(current reference particle)--*/
        	for (j=NCHN; j<N; j++) {	/*--... and accumulate sums--*/
                	if (j==i) continue;
                	rx = fabs(xi-x[j]);
                	ry = fabs(yi-y[j]);
                	rz = fabs(zi-z[j]);
                	while (rx>Lx/2.0) rx -= Lx;
                	while (ry>Ly/2.0) ry -= Ly;
                	while (rz>Lz/2.0) rz -= Lz;
                	r2 = rx*rx + ry*ry + rz*rz;
                	if (r2<rgmax2) {
                        	r = sqrt(r2);
                        	rindx = (int)((r-1.0)/DR);
				if (rindx<NG) n_sum[nblk][rindx]+=1.0;
                	}
		}
 	}  
}
 
void gr_output(vp)	/*--Construct and output final results--*/
double vp;
{
	int i, nblk;
	double g[NBLK][NG], gcont, g_avg[NG], g_var[NG];
	double p_avg, p_var, press[NBLK], diff, p_CS, r;

	/*--zero values--*/
        for (i=0; i<NG; i++) {
		 g_avg[i]=0.0; g_var[i]=0.0;
	}
	p_avg=0.0; p_var=0.0;

	/*--Construct averages for each block--*/
	for (nblk=0; nblk<NBLK; nblk++) {
        	for (i=0; i<NG; i++) {	/*--construct g(r)=navg[i]/n_unc[i]--*/
                	g[nblk][i] = n_sum[nblk][i]/(n_unc[i]*(double)(NSOL*ncnt[nblk]));
			g_avg[i] += g[nblk][i];
		}
        	gcont = 1.5*g[nblk][0]-0.5*g[nblk][1];	/*--linear extrap for g_contact--*/
        	press[nblk] = 1.0 + 4.0*vp*gcont;
		p_avg += press[nblk];
	}
	p_avg /= NBLK;
        for (i=0; i<NG; i++) g_avg[i] /= NBLK;

	/*--Construct variances--*/
	p_var=0.0;
        for (i=0; i<NG; i++) g_var[i]=0.0;
	for (nblk=0; nblk<NBLK; nblk++) {
		diff = p_avg-press[nblk];
		p_var += diff*diff;
        	for (i=0; i<NG; i++) {
			diff = g_avg[i]-g[nblk][i];
			g_var[i] += diff*diff;
		}
	}
	p_var/=(NBLK-1);
        for (i=0; i<NG; i++) g_var[i]/=(NBLK-1);

	/*--Carnahan-Starling expression for pressure--*/
	p_CS = (1.0+vp*(1.0+vp*(1.0-vp)))/((1.0-vp)*(1.0-vp)*(1.0-vp));

	printf("\n P/rho= %5.3f +/- %5.3f\n",p_avg, sqrt(p_var));
	printf("(CS Result: P/rho = %5.3f)\n\n",p_CS);
        for (i=0; i<NG; i++) {	/*--r-value for g(r)=center of g[i] bin--*/
		r=1.0+((double)i+0.5)*DR;
		if (r>sqrt(rgmax2)) break;
                printf("%5.4f  %5.3f +/- %5.3f\n", r, g_avg[i], sqrt(g_var[i]));
	}
}
