/*--------------------------------------------------------------*/
/* File:       hchn_nmer_fld_bigmc.c				*/
/* Usage:      hchn_nmer_fld_bigmc vf				*/
/*                                                              */
/* Monte Carlo Simulation of a Hard-Sphere Flexible n-mer Fluid.*/
/* Algorithm used is intended for small n (now tested for n<=8).*/
/* Program constructs average intermolecular site-site g(r).	*/
/* (Validated for n = 2, 4, 8 [mainly Yethiraj and Hall data]).	*/
/*								*/
/* MC move set: Rigid translation of entire chain.		*/
/*              End-bead rotation about pentultimate bead.	*/
/*		Axial rotation of single interior bead.		*/
/*								*/
/* System consists of NMOL hard-sphere n-mers with bond length	*/
/* L in a rectangular box of volume Vbox=(Lx)(Ly)(Lz) with	*/
/* periodic boundary conditions (pbc).				*/
/* The sphere diameter sigma = unit unit of length.		*/
/* vf = volume fraction of particles;				*/
/* rho = NMOL/Vbox = chain number density;			*/
/* N = n*NMOL = number of spheres in the system;		*/
/* LEN = n = chain length (LEN >= 2).				*/
/*								*/
/* Starting configuration is regular array of x-oriented n-mers.*/
/* System is equilibrated for NEQL MC cycles (1 MC cycle = N	*/
/* attempted moves) followed by NRUN cycles of data production.	*/
/* The max MC displacements for translation (delta), end-bead 	*/
/* rotation (gamma) and axial rotation (alpha) are adjusted	*/
/* during equilibration to give acceptance fractions of about 	*/
/* 50%. Data is collected in NBLK blocks for final statistics.	*/
/* Constructs the radial distribution function g(r) which is	*/
/* calculated out to a distance of min(1+NG*DR,Lx/2,Ly/2,Lz/2).	*/
/*								*/
/* This "bigmc" version uses a neighbor list which is updated	*/
/* every NUPDATE MC cycles.  For this to be efficient, NUPDATE	*/
/* should be set in range 5-20 and the resulting distance rlist	*/
/* must be less than half the box diagonal. List updating is	*/
/* turned off when rlist exceeds this max size.			*/
/*								*/
/*							6/10	*/
/*--------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DIAG	FALSE	/*-- Flag for diagnostic output --*/
#define GFLG	FALSE

#define NCHN	3 	/*-- Chain Length of the Chain--*/
#define	LEN	2	/*-- Solvent-Chain length (LEN >= 2) --*/
#define L	1.0	/*-- Chain bond length (L<=1.0) --*/
#define NMOL	128    	/*-- Total number of n-mers in system --*/

#define NRUN	1000	/*-- # of production MC cycles --*/

#define NG      (200)	/*-- Number of bins for g(r) --*/
#define DR	0.020   /*-- g(r) bin width --*/
#define NBIN	20	/*--number of bins per sigma for w1N(r)--*/
#define NW	( (NCHN-2)*NBIN )  /*-total number of bins for w1N(r)-*/


/*-----------------USUALLY NO CHANGES BELOW HERE-------------------*/
#define	NBND	(LEN-1)		/*--Number of bonds in chain--*/
#define N	(LEN*NMOL)      /*-- Total number of solvent spheres in system --*/
#define NUM	N+NCHN		/*-- Total number of spheres in the system --*/

#define NEQL    (2*NRUN/NBLK)   /*-- # of equil MC cycles = 2 blocks --*/
#define NBLK	10	/*-- # of blocks for statistics --*/

#define NUPDATE	10	/*-- Update time for neighbor list --*/

#define NACCUM	10	/*-- # of MC cycles between data accumulation --*/
#define NADJUST	20	/*-- # of attempts to adjust delta during equil --*/

/*--put rlist def here to make sure it is same everywhere--*/
#define RLIST(delta)	(1.0 + 2.0*NUPDATE*delta + (LEN/2)*L)

#ifndef ISEED		/*-- seed for random num generator --*/
#define ISEED	12761	/*--usually defined on compile line in Makefile--*/
#endif

#define TRUE    1
#define FALSE   0

/*-- uniform random # in range -1.0->1.0 --*/
#define NORM    2147483647.0    /*--(2^31)-1 = largest integer--*/
#define RND()   ( (1.0-2.0*((double)random()/NORM)) )

/*-----Global Variables-----*/
int ncnt[NBLK];		/*--number of data points taken in each block--*/
int lst_ptr[N], nlist[N*N]; /*--neighbor list and pointer arrays--*/
double x[NUM], y[NUM], z[NUM]; /*--particle coordinates--*/
double r2_sum[NBLK],Rg_sum[NBLK],n_sum[NBLK][NG],w1N_sum[NBLK][NW];
double	n_sum[NBLK][NG],/*--g(r) accumulator for each block--*/
	n_unc[NG],	/*--g(r) histogram for ideal gas (uncorrelated) system--*/
	rgmax2;		/*--square of max distance for g(r)--*/
double Lx, Ly, Lz;	/*--simulation box dimensions--*/
double L2, rlist2;	/*--square of bond length and n-list radius--*/
double pi=M_PI;
int overlap();		/*--routine that checks for hard-core overlaps--*/
void chk_bnds();	/*--routine that checks bondlengths--*/
double SIGMA;
int chn_overlap();
void chkchn_bnds();


main(argc,argv)
int argc; char *argv[];
{
        int i, j, icnt, max_list_size, ncyc, nblk;
	int cx, cy, cz, M, ci, isite,k,choice;/* chain part addition */
        int nx, ny, nz, nx_max, ny_max, nz_max;
	int imove, site1, site2, site_flg, chn_num, update_cnt;
	int translate(), end_rot(), axl_rot();
	int chn_trans(),chn_end(),chn_axl();
	double g0, g1, gcont;
        double vf, vcp, rho, vcap, vchn, Vbox, Lbox, wx, wy;
        double delta, gamma, alpha, delta1, alpha1, gamma1;
        double trn_cnt, trn_accpt, f_trn_accpt;
	double chaintrans_cnt, chainend_cnt, chainaxl_cnt;
	double chntrans_accpt, chnend_accpt, chnaxl_accpt;/*Chain part addition*/
        double rot_cnt, rot_accpt, f_rot_accpt;
        double axl_cnt, axl_accpt, f_axl_accpt;
	double chntrans,chnaxl,chnend;
        double rlist, rlist_max,dr,xdel,ydel,zdel;
	double r, r2, r2_min, r2_max, r2_avg,Rg_avg;
	double w1N[NBLK][NW],w1N_avg[NW];
	void list_update(), dump_config(), g_accum(), gr_setup(), gr_output();
	
	/*--get volume fraction from command line--*/
	if (argc!=2) {fprintf(stderr,"\nUsage: %s vf\n\n", argv[0]); exit(1);}
	vf = fabs( atof(argv[1]) );

	vcp=pi/(3.0*sqrt(2.0));	/*--close packed volume fraction (for hard spheres)--*/
	if (vf>vcp) {fprintf(stderr,"\nvf too large\n\n"); exit(1);}
	if (vf>0.48) printf("\nWarning: vf exceeds hard-sphere fluid phase stability limit\n");

	SIGMA=L;
	L2=L*L;
	dr=SIGMA/(double)NBIN;
	
	/*--compute molecular volume and number density--*/
	vcap = (pi/12.0)*(1.0-1.5*L+0.5*L*L2);	/*--volume of spherical cap of height (sigma-L)/2--*/
	vchn = LEN*pi/6.0 - 2.0*(LEN-1)*vcap;	/*--volume of LEN-mer chain with bondlength L--*/
	rho = vf/vchn;	/*--particle number density = vf/Vparticle--*/

	Vbox=(double)N/rho;	/*--box volume--*/
	Lbox=cbrt(Vbox);	/*--edge length for cubic box--*/

        delta=cbrt(vcp/vf)-1.0;	/*--initial max MC x,y,z displacement--*/
	gamma = gamma1=0.5;		/*--initial max MC rotation param--*/
	alpha = alpha1=0.5;		/*--initial max MC axl_rot param--*/
	delta1=delta;

	max_list_size = NMOL*(LEN*LEN*(NMOL-1) + 2*(LEN-2) + (LEN-2)*(LEN-3));
	rlist=RLIST(delta);	/*--neighbor list radius--*/
	rlist2 = rlist*rlist;	/*--square of n-list radius--*/

	/*--Initial Output of Run Details--*/
	printf("\nMonte Carlo Results for L=%4.2f Hard-Sphere %d-mer Fluid with N=%d at vf=%6.5f\n", 
											L, LEN, NMOL, vf);
	printf("Neighbor list in use with NUPDATE=%d\n", NUPDATE);
        printf("Lbox/sigma=%4.2f (rho=%6.5f) DR=%4.3f\n", Lbox, rho, DR);
	printf("ISEED=%d NEQL=%4.2e NRUN=%4.2e NACCUM=%d\n",
					ISEED,(double)NEQL,(double)NRUN,NACCUM);
        printf("\ninitial delta=%4.2f  initial rlist=%3.2f\n", delta, rlist);
	fflush(stdout);

/*----- Setup Initial Config: layers of x-oriented LEN-mers (x,y,z)-spacing = wx+NBND*L,wy,wy --*/
	wx = 1.01;	/*--this value allows for pretty high density--*/
	wy = 2.00;	/*--this value gets down-adjusted below--*/
        nx_max = (int)(Lbox/(wx+NBND*L)+0.5);	/*--num of chns in x-direc--*/
        ny_max = (int)(Lbox/wy);		/*--num of chns in y-direc--*/
	Lx=(double)nx_max*(wx+NBND*L);		/*--starting box dimensions--*/
	Ly=(double)ny_max*(wy); Lz=Vbox/(Lx*Ly);
        nz_max = (int)(Lz/wy);			/*--num of chns in z-direc--*/
	while (nx_max*ny_max*nz_max<NMOL) { 	/*--adjust wy so all chains will fit--*/
		wy = 1.0 + (wy-1.0)/2.0;	/*--reduce wy--*/
        	ny_max = (int)(Lbox/wy);
		Ly=(double)ny_max*(wy); Lz=Vbox/(Lx*Ly);  /*--adjust Ly, Lz--*/
        	nz_max = (int)(Lz/wy);
		if (wy<1.01) {	/*--quit if spacing is too small--*/
			printf("Unable to adjust box ... try changing particle number to: %d\n",
										nx_max*ny_max*nz_max);
			exit(1);
		}
	}
	printf("Lx=%4.2f Ly=%4.2f Lz=%4.2f (nx_max=%d ny_max=%d nz_max=%d)\n",
						 Lx,Ly,Lz,nx_max,ny_max,nz_max); fflush(stdout);
/*----placing the array of the nmer-fluid particles.----*/
        for (i=0, nz=0; i<N; nz++) {	/*--build regular array (rows/layers are not staggered)--*/
	    for (ny=0; (ny<ny_max && i<N); ny++) {
		for (nx=0; (nx<nx_max && i<N); nx++) { /*-insert line of LEN-mers-*/
			x[i]=nx*(wx+NBND*L); y[i]=ny*wy; z[i]=nz*wy;	/*--1st sphere--*/
			for (j=1; j<LEN; j++) {
				x[i+j]=x[i+j-1]+L; y[i+j]=y[i]; z[i+j]=z[i];	/*--(j+1)th sphere--*/
			}
			i+=LEN;	/*--increment sphere counter--*/
		}
	    }
	}

/*--- now the chain part ------------*/
	ci=i;cz=nz-1;cx=nx-1;cy=ny-1;

	if(nz<nz_max){	/* for the new line of in z direction */
	      	x[i]=0.0; y[i]=0.0; z[i]=nz*wy;	
		for (j=1; j<NCHN; j++) {
			x[i+j]=x[i+j-1]+L; y[i+j]=y[i]; z[i+j]=z[i];	/*--(j+1)th sphere--*/
		}	
	}else{
		if (ny<ny_max){
			x[i]=0.0; y[i]=ny*wy;; z[i]=cz*wy;	
			for (j=1; j<NCHN; j++) {
			x[i+j]=x[i+j-1]+L; y[i+j]=y[i]; z[i+j]=z[i];	/*--(j+1)th sphere--*/
			}	
		}else{	
			if( (Lx-(cx*(wx+NBND*L)))>=((NCHN-1)*L+0.5)){
		 		x[i]=cx*(wx+NBND*L); y[i]=(cy)*wy; z[i]=cz*wy;	/*--1st sphere--*/
				for (j=1; j<NCHN; j++) {
					x[i+j]=x[i+j-1]+L; y[i+j]=y[i]; z[i+j]=z[i];	/*--(j+1)th sphere--*/
				}
			}else{
				printf("UNABLE TO ADJUST THE PARTICLES IN THE BOX!!!");
				exit(1);
			}
	    	}
	}
/*------ chain placed ----------*/
	for( i=N;i<NUM;i++){	/*checking Chain overlap*/
		if(chn_overlap(i,0,NUM)==TRUE){printf("there is an overlap in the chain setup\n");fflush(stdout);	
		exit(1);
		}
	}

	if (DIAG) dump_config(0);

	chk_bnds();	/*--verify bond-lengths--*/
	chkchn_bnds();

	rlist_max = sqrt(Lx*Lx+Ly*Ly+Lz*Lz)/2.0;	/*--rlist_max = half box diagonal--*/
	list_update();		/*--Initial setup of neighbor list--*/
	
	for (i=0; i<N; i++) {	/*-- Check initial config for the solvent particles --*/
                if (overlap(i) == TRUE) {
                        printf("Overlapping initial configuration!\n");
			dump_config(0);
                        exit(1);
                }
	}

	if (DIAG) dump_config(0);
	if (GFLG)gr_setup(rho);	/*--Initializations for g(r) calc--*/
        
	srandom(ISEED);	/*--initialize random num generator--*/
/*----------- SYSTEM EQULIBRATION -----------*/
        
	trn_accpt=0.0;rot_accpt=0.0;axl_accpt=0.0;
        trn_cnt=0.0;rot_cnt=0.0;axl_cnt=0.0;
	update_cnt=0;
	chaintrans_cnt=0.0;chainend_cnt=0.0;chainaxl_cnt=0.0;
	chntrans_accpt=0.0;chnend_accpt=0.0;chnaxl_accpt=0.0;
	f_axl_accpt=0.0;	
	chntrans=chnend=chnaxl=0.0;

        for (ncyc=0; ncyc<NEQL; ncyc++) {	/*--perform NEQL MC cycles--*/
		if ( rlist<rlist_max && (ncyc+1)%NUPDATE==0 ) list_update();
		
/*----------- 1 MC cycle = N=LEN*NMOL+NCHN attempted solvent moves-------------*/
                for (icnt=0; icnt<NUM; icnt++) {
		
			choice = random()%NUM;
			if (choice<N){
/*-------------------Fluid Moves--------------------*/
				imove = random()%2;	/*--select move type (0=trans, 1=rot)--*/
				site1 = choice; 	/*--select chain site1 at random--*/
				chn_num = site1/LEN;	/*--parent chain (note use of integer truncation)--*/
				site_flg = site1%LEN;	/*--site position along chain--*/

				if (imove==1) {	/*--setup for the rotation moves--*/
					switch(site_flg) {
						case 0:     site2=site1+1; break; /*--site1 = left end--*/
						case LEN-1: site2=site1-1; break; /*--site1 = right end--*/
						default:    imove=2; break;  /*--axial rotation for middle sites--*/
					}
				}
				switch(imove) {
					case 0: trn_cnt++; /*--Translation Move Attempt--*/
						if (translate(chn_num,delta)==TRUE)
							trn_accpt++; /*--increment accept counter--*/
						break;	
					case 1: rot_cnt++; /*--Endbead-Rot Attempt (site1=end, site2=connected)--*/
						if (end_rot(site1,site2,gamma)==TRUE)
							rot_accpt++;    /*--increment accept counter--*/
						break;
					case 2: axl_cnt++;	/*--Axial-Rot Attempt (site1=non-end)--*/
						if (axl_rot(site1,alpha)==TRUE)
							axl_accpt++;    /*--increment accept counter--*/
						break;
					default: printf("Error: Bad imove value: %d\n", imove); fflush(stdout);
				}
			chk_bnds();
/*---------------------- Fluid Moved-------------------*/
	        	}else{
/*-----------------------Chain moves----------------*/
				imove = random()%2;	/*--select move type (0=trans, 1=rot)--*/
				isite = choice;		/*-- selected the isite position of the chain--*/	

				if (imove==1) {	/*--setup for the rotation moves--*/
					switch(isite) {
						case N:     site2=isite+1; break; /*--isite= left end--*/
						case NUM-1: site2=isite-1; break; /*--isite = right end--*/
						default:    imove=2; break;  /*--axial rotation for middle sites--*/
					}
				}
				switch(imove) {
					case 0: chaintrans_cnt++; /*--Translation Move Attempt--*/
					if (chn_trans(delta1)==TRUE){
							chntrans_accpt++; /*--increment accept counter--*/
						}
						break;
					case 1: chainend_cnt++; /*--Endbead-Rot Attempt (site1=end, site2=connected)--*/
						if (chn_end(isite,site2,gamma1)==TRUE){
							chnend_accpt++;    /*--increment accept counter--*/
						}
						break;
					case 2: chainaxl_cnt++;	/*--Axial-Rot Attempt (site1=non-end)--*/
						if (chn_axl(isite,alpha1)==TRUE){
							chnaxl_accpt++;    /*--increment accept counter--*/
						}
						break;
					default: printf("Error: Bad imove value: %d\n", imove); fflush(stdout);
				}
			chkchn_bnds();
			}
		}
/*-------------------Chain moved-------------------*/
	
		if ((ncyc+1)%(NEQL/NADJUST)==0) {	/*--adjust delta, gamma, and alpha--*/
			update_cnt++;		   /*--target = 50% acceptance of each move type--*/

			/*--delta adjust ... for translation--*/
			f_trn_accpt = trn_accpt/trn_cnt;
			if (f_trn_accpt<0.5) delta/=1.1;
			else delta*=1.1;

			rlist=RLIST(delta);  /*--neighbor list radius--*/
			rlist2=rlist*rlist;
			list_update();	/*--force list update on new delta--*/

			/*--delta1 adjust for the chain--*/
			chntrans = chntrans_accpt/chaintrans_cnt;
			if (chntrans<0.5) delta1/=1.1;
			else delta1*=1.1;

			/*--gamma adjust ... for end rotation--*/
			f_rot_accpt = rot_accpt/rot_cnt;
			if (f_rot_accpt<0.5) gamma/=1.1;
			else gamma*=1.1;
			if (gamma>=1.0) gamma=1.0;	/*--upper limit on gamma (best choice?)--*/

			/*-----gamma1 adjust for the chain part-------*/
			chnend = chnend_accpt/chainend_cnt;
			if (chnend<0.5) gamma1/=1.1;
			else gamma1*=1.1;
			if (gamma1>=1.0) gamma1=1.0;	/*--upper limit on gamma (best choice?)--*/

			
			/*--alpha adjust .. for interior bead axial rotation--*/
			if (LEN>2.0) {
				f_axl_accpt = axl_accpt/axl_cnt;
				if (f_axl_accpt<0.5) alpha/=1.1;
				else alpha*=1.1;
				if (alpha>pi){ alpha=pi;}	/*--upper limit on alpha--*/
			}

			if ((NCHN>2) && (chainaxl_cnt>0)){/*-------alpha1 adjust for the chain part-----*/
				chnaxl=chnaxl_accpt/chainaxl_cnt;
				if (chnaxl<0.5) alpha1/=2.1;
				else alpha1*=1.1;
				if (alpha1>pi) alpha1=pi;
			}

			printf("Equil Update %d: delta=%5.4f (rlist=%5.4f) gamma=%5.4f alpha=%5.4f ",
								update_cnt, delta, rlist, gamma, alpha);
			printf(" nlist_end=%d (max=%d)\n", lst_ptr[N], max_list_size);
			printf("            : delta1=%5.4f gamma1=%5.4f alpha1=%5.4f \n",
								delta1, gamma1, alpha1);
			fflush(stdout);
        
			/*--re-zero counters for next equil step--*/
			trn_cnt=rot_cnt=axl_cnt=0.0;
			trn_accpt=rot_accpt=axl_accpt=0.0;
			chaintrans_cnt=chainend_cnt=chainaxl_cnt=0.0;
			chntrans_accpt=chnend_accpt=chnaxl_accpt=0.0;
			
			chk_bnds();
			chkchn_bnds();  /*----verify bond lengths of the chain through equilibration----*/
		}
         }
        printf("NEQL= %3.2e  (Final Equil tran_accept: ", (double)NEQL);
	printf(" %3.1f%%  rot_accept: %3.1f%%  axl_accept: %3.1f%%)\n",
					 100.0*f_trn_accpt, 100.0*f_rot_accpt, 100.0*f_axl_accpt);
        printf("               (Final Equil chain part  chntrans_accpt: " );
	printf(" %3.1f%%  chnend_accpt: %3.1f%%  chnaxl_accpt: %3.1f%%)\n\n",
					 100.0*chntrans, 100.0*chnend, 100.0*chnaxl);
	if (rlist > rlist_max) printf("rlist exceeds rlist_max (=%5.3f): No list updating\n", rlist_max);
	fflush(stdout);

	if (DIAG) dump_config(NEQL);
/*------------- DATA PRODUCTION -----------------*/
	trn_accpt=0.0;rot_accpt=0.0;axl_accpt=0.0;
        trn_cnt=0.0;rot_cnt=0.0;axl_cnt=0.0;
	update_cnt=0;
	chaintrans_cnt=0.0;chainend_cnt=0.0;chainaxl_cnt=0.0;
	chntrans_accpt=0.0;chnend_accpt=0.0;chnaxl_accpt=0.0;
	f_axl_accpt=0.0;	
	chntrans=chnend=chnaxl=0.0;

	for (ncyc=0; ncyc<NRUN; ncyc++) {	/*--perform NRUN MC cycles--*/
		if ( rlist<rlist_max && (ncyc+1)%NUPDATE==0 ) list_update();
		
/*----------- 1 MC cycle = N=LEN*NMOL+NCHN attempted solvent moves-------------*/

/*----------- 1 MC cycle = N=LEN*NMOL+NCHN attempted solvent moves-------------*/
                for (icnt=0; icnt<NUM; icnt++) {
		
			choice = random()%NUM;
			if (choice<N){
/*-------------------Fluid Moves--------------------*/
				imove = random()%2;	/*--select move type (0=trans, 1=rot)--*/
				site1 = choice; 	/*--select chain site1 at random--*/
				chn_num = site1/LEN;	/*--parent chain (note use of integer truncation)--*/
				site_flg = site1%LEN;	/*--site position along chain--*/

				if (imove==1) {	/*--setup for the rotation moves--*/
					switch(site_flg) {
						case 0:     site2=site1+1; break; /*--site1 = left end--*/
						case LEN-1: site2=site1-1; break; /*--site1 = right end--*/
						default:    imove=2; break;  /*--axial rotation for middle sites--*/
					}
				}
				switch(imove) {
					case 0: trn_cnt++; /*--Translation Move Attempt--*/
						if (translate(chn_num,delta)==TRUE)
							trn_accpt++; /*--increment accept counter--*/
						break;	
					case 1: rot_cnt++; /*--Endbead-Rot Attempt (site1=end, site2=connected)--*/
						if (end_rot(site1,site2,gamma)==TRUE)
							rot_accpt++;    /*--increment accept counter--*/
						break;
					case 2: axl_cnt++;	/*--Axial-Rot Attempt (site1=non-end)--*/
						if (axl_rot(site1,alpha)==TRUE)
							axl_accpt++;    /*--increment accept counter--*/
						break;
					default: printf("Error: Bad imove value: %d\n", imove); fflush(stdout);
				}
printf("imove=%i\n",imove);fflush(stdout);
			chk_bnds();
/*---------------------- Fluid Moved-------------------*/
	        	}else{
/*-----------------------Chain moves----------------*/
				imove = random()%2;	/*--select move type (0=trans, 1=rot)--*/
				isite = choice;		/*-- selected the isite position of the chain--*/	

				if (imove==1) {	/*--setup for the rotation moves--*/
					switch(isite) {
						case N:     site2=isite+1; break; /*--isite= left end--*/
						case NUM-1: site2=isite-1; break; /*--isite = right end--*/
						default:    imove=2; break;  /*--axial rotation for middle sites--*/
					}
				}
				switch(imove) {
					case 0: chaintrans_cnt++; /*--Translation Move Attempt--*/
					if (chn_trans(delta1)==TRUE){
							chntrans_accpt++; /*--increment accept counter--*/
						}
						break;
					case 1: chainend_cnt++; /*--Endbead-Rot Attempt (site1=end, site2=connected)--*/
						if (chn_end(isite,site2,gamma1)==TRUE){
							chnend_accpt++;    /*--increment accept counter--*/
						}
						break;
					case 2: chainaxl_cnt++;	/*--Axial-Rot Attempt (site1=non-end)--*/
						if (chn_axl(isite,alpha1)==TRUE){
							chnaxl_accpt++;    /*--increment accept counter--*/
						}
						break;
					default: printf("Error: Bad imove value: %d\n", imove); fflush(stdout);
				}
			chkchn_bnds();
			}
		}
/*-------------------Chain moved-------------------*/



/* Now the calculation after each accumulation cycle of each block */
		if (ncyc%NACCUM==0) {      /*--take data every NACCUM cycle-*/
			nblk=ncyc/(NRUN/NBLK); 	/*--current block number--*/
			if (nblk>=NBLK) {printf("Error: nblk too large!\n"); fflush(stdout); }
			ncnt[nblk]++;		/*--increment counter--*/
			if (GFLG) g_accum(nblk);		/*--accumulate data--*/
		
/*			xdel=x[NCHN+N-1]-x[0+N]; ydel=y[NCHN+N-1]-y[0+N]; zdel=z[NCHN+N-1]-z[0+N];
                        r2 = (xdel*xdel+ydel*ydel+zdel*zdel);
			if (r2<r2_min) r2_min=r2;
			if (r2>r2_max) r2_max=r2;
                        r2_sum[nblk] += r2;  */   /*--end-to-end distance--*/
			
			/*--build histogram for P1N(r)--*/
			k=(int)((sqrt(r2)-SIGMA)/dr);
			if (k<NW) w1N_sum[nblk][k]++;

			for (i=N; i<NCHN+N-1; i++){
                            for (j=i+1; j<NCHN+N; j++) {/*-radius of gyration-*/
                                      xdel=x[j]-x[i]; ydel=y[j]-y[i]; zdel=z[j]-z[i];
                                      r2 = (xdel*xdel+ydel*ydel+zdel*zdel);
                                      Rg_sum[nblk] += r2/(double)(NCHN*NCHN);
                     	}
			}

        	       /*--Output data after each block complete--*/
   	               if (ncnt[nblk]==(NRUN/(NBLK*NACCUM))) {
   				if (GFLG){
		             		g0 = n_sum[nblk][0]/(n_unc[0]*(double)(N*ncnt[nblk]));
	                		g1 = n_sum[nblk][1]/(n_unc[1]*(double)(N*ncnt[nblk]));
	        			gcont = 1.5*g0-0.5*g1;	/*--linear extrap for g_contact--*/
	                                printf("Block %d complete (%d cycles): ", nblk+1, (nblk+1)*NRUN/NBLK);
	                                printf("g_contact= %6.4f\n", gcont);
	                                fflush(stdout);
				}else{
					printf("       Block %d complete (%d cycles): \n", nblk+1, (nblk+1)*NRUN/NBLK);
					fflush(stdout);
        			}
          		}	
		
		chk_bnds(); /*--verify bond lengths after each block--*/
		chkchn_bnds();/*--verify bond lengths of the chain after each block completed--*/		
		}	
	}
/*----Data production complete----------*/
	for( i=N;i<NUM;i++){	
		if(chn_overlap(i,0,NUM)==TRUE){printf("there is an overlap in the chain \n");fflush(stdout);	
		exit(1);
		}
	}

	f_trn_accpt = trn_accpt/trn_cnt;
	chntrans = chntrans_accpt/chaintrans_cnt;
	f_rot_accpt = rot_accpt/rot_cnt;
	chnend = chnend_accpt/chainend_cnt;
	printf("NRUN= %3.2e  (Final Prod accept for the solvent: ", (double)NRUN);
       printf("Translate=%3.1f%%, End-bead rotation=%3.1f%%", 
				100.0*f_trn_accpt, 100.0*f_rot_accpt);

	if (LEN>2.0) {
		f_axl_accpt = axl_accpt/axl_cnt;
		printf( "Axial Rotation=%3.2f%%",100.0*f_axl_accpt);
	}
	printf("\nNRUN= %3.2e  (Final Prod accept for the Chain: ", (double)NRUN);
        printf("Translate=%3.1f%%, End-bead rotation=%3.1f%%", 
				100.0*chntrans, 100.0*chnend);

	if ((NCHN>2) && (chainaxl_cnt>0)){
		chnaxl=chnaxl_accpt/chainaxl_cnt;
		printf( "Axial Rotation=%3.2f%%",100.0*chnaxl);
	}

/*-----Data output and results-------------*/
	if (GFLG) gr_output(vf);
	
	r2_avg=Rg_avg=0.0;
        for (i=0; i<NW; i++) {
		w1N_avg[i]=0.0;
	}

	for (nblk=0; nblk<NBLK; nblk++) {
              r2_avg += r2_sum[nblk]/ncnt[nblk];
              Rg_avg += Rg_sum[nblk]/ncnt[nblk];
        	for (i=0; i<NW; i++) {	/*--construct P1N(r)--*/
                	w1N[nblk][i] = (SIGMA/dr)*w1N_sum[nblk][i]/(ncnt[nblk]);
			w1N_avg[i] += w1N[nblk][i];
		}
        }





}

int chn_trans(delta1) /*--rigid translation of single chain--*/
double delta1;
{
	
	int i, site,  ovrlap_flg; 
	double cxsave[NCHN], cysave[NCHN], czsave[NCHN];
	double dx, dy, dz;
	void pbc_adjust();

/*--save site positions--*/
	for (i=0; i<NCHN; i++) {
		site = i+N;
		cxsave[i]=x[site]; cysave[i]=y[site]; czsave[i]=z[site];
	}

	ovrlap_flg=FALSE;	/*--initialize flag for no overlap--*/
	dx = delta1*RND(); dy = delta1*RND(); dz = delta1*RND();	/*--random displacement vector--*/
	for (i=0; i<NCHN; i++) {	/*--translate all beads ... --*/
		site = i+N;
		x[site]+=dx; y[site]+=dy; z[site]+=dz;
	}
	for (i=0; i<NCHN; i++) {	/*-- ... then check for overlap--*/
		site = i+N;
		pbc_adjust(site);
		if (chn_overlap(site,0,NUM)==TRUE) {/*--break out of loop on overlap--*/
			ovrlap_flg=TRUE;
			break;
    		}
	}
	if (ovrlap_flg==TRUE) {
		for (i=0; i<NCHN; i++) {	 	/*--restore positions and  --*/
			site = i+N;
			x[site]=cxsave[i]; y[site]=cysave[i]; z[site]=czsave[i];
		}
		return(FALSE);			/*-- ...  return FALSE on overlap--*/
	} else {
		for (i=0; i<NCHN; i++) {		/*--cleanup and  --*/
			site = i+N;
			pbc_adjust(site);
		}
		return(TRUE);			/*-- ... return TRUE on success--*/
	}

}

int translate(chn_num,delta) /*--rigid translation of single chain--*/
int chn_num; double delta;
{
	int i, site, ref_site, ovrlap_flg;
	double xsave[LEN], ysave[LEN], zsave[LEN];
	double dx, dy, dz;
	void pbc_adjust();

	ref_site = chn_num*LEN;  /*--index of 1st site in chain chn_num--*/

	/*--save site positions--*/
	for (i=0; i<LEN; i++) {
		site = i+ref_site;
		xsave[i]=x[site]; ysave[i]=y[site]; zsave[i]=z[site];
	}

	ovrlap_flg=FALSE;	/*--initialize flag for no overlap--*/
	dx = delta*RND(); dy = delta*RND(); dz = delta*RND();	/*--random displacement vector--*/
	for (i=0; i<LEN; i++) {	/*--translate all beads ... --*/
		site = i+ref_site;
		x[site]+=dx; y[site]+=dy; z[site]+=dz;
	}
	for (i=0; i<LEN; i++) {	/*-- ... then check for overlap--*/
		site = i+ref_site;
		if (overlap(site)==TRUE) {/*--break out of loop on overlap--*/
			ovrlap_flg=TRUE;
			break;
		}
	}
	if(ovrlap_flg==FALSE){
		for( i=N;i<NUM;i++){	/*checking Chain overlap*/
			if(chn_overlap(i,0,NUM)==TRUE){
			ovrlap_flg=TRUE;	
			break;}
		}
	}

	if (ovrlap_flg==TRUE) {
		for (i=0; i<LEN; i++) {	 	/*--restore positions and  --*/
			site = i+ref_site;
			x[site]=xsave[i]; y[site]=ysave[i]; z[site]=zsave[i];
		}
		return(FALSE);			/*-- ...  return FALSE on overlap--*/
	} else {
		for (i=0; i<LEN; i++) {		/*--cleanup and  --*/
			site = i+ref_site;
			pbc_adjust(site);
		}
		return(TRUE);			/*-- ... return TRUE on success--*/
	}
}

int end_rot(rotsite,fixsite,gamma) /*--rotates rotsite(=end bead) about fixsite--*/
int rotsite, fixsite; double gamma; /*--See Frenkel and Smit, sec 3.3.2--*/
{
	int i,ovrlap_flg;
	double ux, uy, uz, u, u2;
	double vx, vy, vz, vxnew, vynew, vznew, vnew;
	double x1, y1, z1, x2, y2, z2, xdel, ydel, zdel;
	void pbc_adjust();

	/*--save site positions--*/
	x1=x[fixsite]; y1=y[fixsite]; z1=z[fixsite];
	x2=x[rotsite]; y2=y[rotsite]; z2=z[rotsite];

	/*--make coordinate adjustments for chain "broken" by pbc--*/
        while ((xdel=fabs(x1-x2))>Lx/2.0) {if (x1>x2) x1 -= Lx; else x1 += Lx;}
        while ((ydel=fabs(y1-y2))>Ly/2.0) {if (y1>y2) y1 -= Ly; else y1 += Ly;}
        while ((zdel=fabs(z1-z2))>Lz/2.0) {if (z1>z2) z1 -= Lz; else z1 += Lz;}

	/*-generate random vector on unit sphere-*/
	do { ux=RND(); uy=RND(); uz=RND();
	} while ((u2=ux*ux+uy*uy+uz*uz)>1.0);
	u=sqrt(u2);

	vx = (x2-x1)/L; vy = (y2-y1)/L; vz = (z2-z1)/L; /*--"dimer" orientation vector fix->rot--*/

	/*--construct new unnormalized orientation vector: vnew = v + gamma*u--*/
	vxnew = vx+gamma*ux/u; vynew = vy+gamma*uy/u; vznew = vz+gamma*uz/u;
	vnew = sqrt(vxnew*vxnew + vynew*vynew + vznew*vznew);	/*--vector magnitude--*/

	/*--coordinates for rotsite after "rotation"--*/
	x[rotsite] = x1 + L*vxnew/vnew;
	y[rotsite] = y1 + L*vynew/vnew;
	z[rotsite] = z1 + L*vznew/vnew;

	for( i=N;i<NUM;i++){	/*checking Chain overlap*/
		if(chn_overlap(i,0,NUM)==TRUE){
		ovrlap_flg=TRUE;	
		break;}
	}

	if ((overlap(rotsite)==TRUE)||(ovrlap_flg==TRUE)) {
		x[rotsite]=x2; y[rotsite]=y2; z[rotsite]=z2;	 /*--restore position and--*/
		return(FALSE);				/*-- ...  return FALSE on overlap--*/
	} else {
		pbc_adjust(rotsite);			/*--cleanup and--*/
 		return(TRUE);				/*-- ... return TRUE on success--*/
	}
}

int chn_end(rotsite,fixsite,gamma) /*--rotates rotsite(=end bead) about fixsite--*/
int rotsite, fixsite; double gamma; /*--See Frenkel and Smit, sec 3.3.2--*/
{
	double ux, uy, uz, u, u2;
	double vx, vy, vz, vxnew, vynew, vznew, vnew;
	double x1, y1, z1, x2, y2, z2, xdel, ydel, zdel;
	void pbc_adjust();

	/*--save site positions--*/
	x1=x[fixsite]; y1=y[fixsite]; z1=z[fixsite];
	x2=x[rotsite]; y2=y[rotsite]; z2=z[rotsite];

	/*--make coordinate adjustments for chain "broken" by pbc--*/
        while ((xdel=fabs(x1-x2))>Lx/2.0) {if (x1>x2) x1 -= Lx; else x1 += Lx;}
        while ((ydel=fabs(y1-y2))>Ly/2.0) {if (y1>y2) y1 -= Ly; else y1 += Ly;}
        while ((zdel=fabs(z1-z2))>Lz/2.0) {if (z1>z2) z1 -= Lz; else z1 += Lz;}

	/*-generate random vector on unit sphere-*/
	do { ux=RND(); uy=RND(); uz=RND();
	} while ((u2=ux*ux+uy*uy+uz*uz)>1.0);
	u=sqrt(u2);

	vx = (x2-x1)/L; vy = (y2-y1)/L; vz = (z2-z1)/L; /*--"dimer" orientation vector fix->rot--*/

	/*--construct new unnormalized orientation vector: vnew = v + gamma*u--*/
	vxnew = vx+gamma*ux/u; vynew = vy+gamma*uy/u; vznew = vz+gamma*uz/u;
	vnew = sqrt(vxnew*vxnew + vynew*vynew + vznew*vznew);	/*--vector magnitude--*/

	/*--coordinates for rotsite after "rotation"--*/
	x[rotsite] = x1 + L*vxnew/vnew;
	y[rotsite] = y1 + L*vynew/vnew;
	z[rotsite] = z1 + L*vznew/vnew;
	
	if (chn_overlap(rotsite,0,NUM)==TRUE)  {
		x[rotsite]=x2; y[rotsite]=y2; z[rotsite]=z2;	 /*--restore position and--*/
		return(FALSE);				/*-- ...  return FALSE on overlap--*/
	} else {
		pbc_adjust(rotsite);			/*--cleanup and--*/
 		return(TRUE);				/*-- ... return TRUE on success--*/
	}
}


int axl_rot(is,alpha) /*--axial rotation of site is about (is-1)-(is+1) axis--*/
int is; double alpha;	/*--this is a "clean" version from Wang-Landau program--*/
{
	int i, ovrlap_flg;
        double xp1,xp2,yp1,yp2,zp1,zp2,xi,yi,zi;
        double xdel,ydel,zdel,r,cos_a,sin_a;
	double angle, xsav, ysav, zsav;
	double c00,c01,c02,c10,c11,c12,c20,c21,c22;
	double xfct,yfct,zfct;
	void pbc_adjust();

	/*--save site position--*/
	xsav=x[is]; ysav=y[is]; zsav=z[is];

	angle = alpha*RND(); /*--select rotation angle: -angle -> angle rads--*/

        cos_a=cos(angle); 
	sin_a=sqrt(1.0-cos_a*cos_a);
	if (angle < 0.0) sin_a *= -1.0;

        /*--Determine (theta,phi) orientation of the rotation axis--*/
        xp1=x[is-1]; yp1=y[is-1]; zp1=z[is-1];
        xp2=x[is+1]; yp2=y[is+1]; zp2=z[is+1];

	/*--make coordinate adjustments for chain "broken" by pbc--*/
        while (fabs(xp1-xsav)>Lx/2.0) {if (xp1>xsav) xp1 -= Lx; else xp1 += Lx;}
        while (fabs(yp1-ysav)>Ly/2.0) {if (yp1>ysav) yp1 -= Ly; else yp1 += Ly;}
        while (fabs(zp1-zsav)>Lz/2.0) {if (zp1>zsav) zp1 -= Lz; else zp1 += Lz;}

        while (fabs(xp2-xsav)>Lx/2.0) {if (xp2>xsav) xp2 -= Lx; else xp2 += Lx;}
        while (fabs(yp2-ysav)>Ly/2.0) {if (yp2>ysav) yp2 -= Ly; else yp2 += Ly;}
        while (fabs(zp2-zsav)>Lz/2.0) {if (zp2>zsav) zp2 -= Lz; else zp2 += Lz;}

        xdel=xp2-xp1; ydel=yp2-yp1; zdel=zp2-zp1;
        r=sqrt(xdel*xdel+ydel*ydel+zdel*zdel);

	if (fabs(r-2.0*L)<1.0e-6) return(FALSE); /*--return FALSE on straight segment config--*/

	/*------B rotation aligns rot-axis with z-axis------*/
	/*----------(theta,phi) rotation Matrix B-----------*/
	/*  b00=cos_p*cos_t; b01=cos_p*sin_t; b02=(-sin_p); */
	/*  b10=(-sin_t);    b11=cos_t;       b12=0.0;      */
	/*  b20=sin_p*cos_t; b21=sin_p*sin_t; b22=cos_p;    */
	/*--------------------------------------------------*/

	/*--A rotation is about z-axis--*/
        /*---alpha rotation Matrix A----*/
	/*	 cos_a  sin_a  0.0	*/
 	/*	-sin_a  cos_a  0.0	*/
 	/*	 0.0     0.0   1.0	*/
	/*------------------------------*/

        /*-------Full rotation Matrix is  C=(B^T)A(B)--------*/

	xfct = xdel/r; yfct = ydel/r; zfct = zdel/r;

	/*--construct cij matrix elements--*/
	c00 = (xfct*xfct)*(1.0-cos_a) + cos_a;
	c11 = (yfct*yfct)*(1.0-cos_a) + cos_a;
	c22 = (zfct*zfct)*(1.0-cos_a) + cos_a;

	c10 = (xfct*yfct)*(1.0-cos_a) + (-zfct)*sin_a;
	c01 = (xfct*yfct)*(1.0-cos_a) + ( zfct)*sin_a;

	c20 = (xfct*zfct)*(1.0-cos_a) + ( yfct)*sin_a;
	c02 = (xfct*zfct)*(1.0-cos_a) + (-yfct)*sin_a;

	c12 = (yfct*zfct)*(1.0-cos_a) + ( xfct)*sin_a;
	c21 = (yfct*zfct)*(1.0-cos_a) + (-xfct)*sin_a;

        /*--rotate site is about ip1-ip2 axis--*/
        xi=xsav-xp1; yi=ysav-yp1; zi=zsav-zp1;
        x[is] = xp1 + c00*xi + c01*yi + c02*zi;
        y[is] = yp1 + c10*xi + c11*yi + c12*zi;
        z[is] = zp1 + c20*xi + c21*yi + c22*zi;
	
	for( i=N;i<NUM;i++){	/*checking Chain overlap*/
		if(chn_overlap(i,0,NUM)==TRUE){
		ovrlap_flg=TRUE;	
		break;}
	}


	if ((overlap(is)==TRUE) || (ovrlap_flg==TRUE)) {
		x[is]=xsav; y[is]=ysav; z[is]=zsav;	 /*--restore position and--*/
		return(FALSE);				/*-- ...  return FALSE on overlap--*/
	} else {
		pbc_adjust(is);			/*--cleanup and--*/
 		return(TRUE);		/*-- ... return TRUE on success--*/
	}
}

int chn_axl(is,alpha) /*--axial rotation of site is about (is-1)-(is+1) axis--*/
int is; double alpha;	/*--this is a "clean" version from Wang-Landau program--*/
{
        double xp1,xp2,yp1,yp2,zp1,zp2,xi,yi,zi;
        double xdel,ydel,zdel,r,cos_a,sin_a;
	double angle, xsav, ysav, zsav;
	double c00,c01,c02,c10,c11,c12,c20,c21,c22;
	double xfct,yfct,zfct;
	void pbc_adjust();

	/*--save site position--*/
	xsav=x[is]; ysav=y[is]; zsav=z[is];

	angle = alpha*RND(); /*--select rotation angle: -angle -> angle rads--*/

        cos_a=cos(angle); 
	sin_a=sqrt(1.0-cos_a*cos_a);
	if (angle < 0.0) sin_a *= -1.0;

        /*--Determine (theta,phi) orientation of the rotation axis--*/
        xp1=x[is-1]; yp1=y[is-1]; zp1=z[is-1];
        xp2=x[is+1]; yp2=y[is+1]; zp2=z[is+1];

	/*--make coordinate adjustments for chain "broken" by pbc--*/
        while (fabs(xp1-xsav)>Lx/2.0) {if (xp1>xsav) xp1 -= Lx; else xp1 += Lx;}
        while (fabs(yp1-ysav)>Ly/2.0) {if (yp1>ysav) yp1 -= Ly; else yp1 += Ly;}
        while (fabs(zp1-zsav)>Lz/2.0) {if (zp1>zsav) zp1 -= Lz; else zp1 += Lz;}

        while (fabs(xp2-xsav)>Lx/2.0) {if (xp2>xsav) xp2 -= Lx; else xp2 += Lx;}
        while (fabs(yp2-ysav)>Ly/2.0) {if (yp2>ysav) yp2 -= Ly; else yp2 += Ly;}
        while (fabs(zp2-zsav)>Lz/2.0) {if (zp2>zsav) zp2 -= Lz; else zp2 += Lz;}

        xdel=xp2-xp1; ydel=yp2-yp1; zdel=zp2-zp1;
        r=sqrt(xdel*xdel+ydel*ydel+zdel*zdel);

	if (fabs(r-2.0*L)<1.0e-6) return(FALSE); /*--return FALSE on straight segment config--*/

	/*------B rotation aligns rot-axis with z-axis------*/
	/*----------(theta,phi) rotation Matrix B-----------*/
	/*  b00=cos_p*cos_t; b01=cos_p*sin_t; b02=(-sin_p); */
	/*  b10=(-sin_t);    b11=cos_t;       b12=0.0;      */
	/*  b20=sin_p*cos_t; b21=sin_p*sin_t; b22=cos_p;    */
	/*--------------------------------------------------*/

	/*--A rotation is about z-axis--*/
        /*---alpha rotation Matrix A----*/
	/*	 cos_a  sin_a  0.0	*/
 	/*	-sin_a  cos_a  0.0	*/
 	/*	 0.0     0.0   1.0	*/
	/*------------------------------*/

        /*-------Full rotation Matrix is  C=(B^T)A(B)--------*/

	xfct = xdel/r; yfct = ydel/r; zfct = zdel/r;

	/*--construct cij matrix elements--*/
	c00 = (xfct*xfct)*(1.0-cos_a) + cos_a;
	c11 = (yfct*yfct)*(1.0-cos_a) + cos_a;
	c22 = (zfct*zfct)*(1.0-cos_a) + cos_a;

	c10 = (xfct*yfct)*(1.0-cos_a) + (-zfct)*sin_a;
	c01 = (xfct*yfct)*(1.0-cos_a) + ( zfct)*sin_a;

	c20 = (xfct*zfct)*(1.0-cos_a) + ( yfct)*sin_a;
	c02 = (xfct*zfct)*(1.0-cos_a) + (-yfct)*sin_a;

	c12 = (yfct*zfct)*(1.0-cos_a) + ( xfct)*sin_a;
	c21 = (yfct*zfct)*(1.0-cos_a) + (-xfct)*sin_a;

        /*--rotate site is about ip1-ip2 axis--*/
        xi=xsav-xp1; yi=ysav-yp1; zi=zsav-zp1;
        x[is] = xp1 + c00*xi + c01*yi + c02*zi;
        y[is] = yp1 + c10*xi + c11*yi + c12*zi;
        z[is] = zp1 + c20*xi + c21*yi + c22*zi;

	if (chn_overlap(is,0,NUM)==TRUE) {
		x[is]=xsav; y[is]=ysav; z[is]=zsav;	 /*--restore position and--*/
		return(FALSE);				/*-- ...  return FALSE on overlap--*/
	} else {
		pbc_adjust(is);			/*--cleanup and--*/
 		return(TRUE);		/*-- ... return TRUE on success--*/
	}
}


void pbc_adjust(i) /*--adjust coordinates of site i in global x,y,z arrays for pbc--*/
int i;
{
	if      (x[i]<0.0) x[i]+=Lx;
        else if (x[i]>Lx ) x[i]-=Lx;
        if      (y[i]<0.0) y[i]+=Ly;
        else if (y[i]>Ly ) y[i]-=Ly;
        if      (z[i]<0.0) z[i]+=Lz;
        else if (z[i]>Lz ) z[i]-=Lz;
}

int chn_overlap(i,nstrt,nstop) /*--check for overlap of particle i--*/
int i,nstrt,nstop;	   /*--PBC's enforced even if chain outside the box--*/
{
        int j;
        double xi,yi,zi,rx,ry,rz,r2,d2,SIGMA2;
	SIGMA2=L*L;
        xi=x[i]; yi=y[i]; zi=z[i];
        for (j=nstrt; j<nstop; j++) {
		
                if (j==i) continue;
		/*--set d2min for current ij pair--*/
		if (i>=N && j>=N) {
			if (fabs(i-j)==1) continue; /*--skip on bonded sites--*/
			else d2=SIGMA2;			/*--chain-chain--*/
		}
		
		else	   d2=SIGMA2;		/*--solvent-chain--*/
                
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

int overlap(i) 	/*--check for overlap of particle i--*/
int i;	   	/*----using neighbor list nlist[]----*/
{
        int j,jl,jlstrt,jlstop;
        double xi,yi,zi,rx,ry,rz,r2;

        xi=x[i]; yi=y[i]; zi=z[i];

	jlstrt=lst_ptr[i]; jlstop=lst_ptr[i+1]; /*--start and stop positions in nlist--*/
        for (jl=jlstrt; jl<jlstop; jl++) {
		j = nlist[jl];	/*--neighbor index--*/
                rx = fabs(xi-x[j]);
                ry = fabs(yi-y[j]);
                rz = fabs(zi-z[j]);
                while (rx>Lx/2.0) rx -= Lx;	/*--pbc adjustments--*/
                while (ry>Ly/2.0) ry -= Ly;
                while (rz>Lz/2.0) rz -= Lz;
                r2 = rx*rx + ry*ry + rz*rz;
                if (r2<1.0) return(TRUE);	/*--overlap here ... stop check--*/
        }
        return(FALSE);	/*--no overlap--*/
}

void list_update()	/*--update neighbor list for LEN-mer chain fluid--*/
{
	int i1, i2, i3, j, ptr_val, site_flg;
	double xi,yi,zi,rx,ry,rz,r2;
	ptr_val=0;
	for (i1=0; i1<N; i1++) {	/*--construct list for site i1--*/
		lst_ptr[i1]=ptr_val;	/*--record start location in list for this site--*/
        	xi=x[i1]; yi=y[i1]; zi=z[i1];
		site_flg = i1%LEN;	/*--set flag for site location on chain--*/
		switch(site_flg) { /*--indices of bonded partners--*/
			case 0:     i2=i1+1; i3=(-1); break;   /*--left end (i3 not used)--*/
			case LEN-1: i2=i1-1; i3=(-1); break;   /*--right end (i3 not used)--*/
			default:    i2=i1-1; i3=i1+1; break; /*--interior: two bonded sites here--*/
		}
		for (j=0; j<N; j++) {	/*--scan through all sites j--*/
			if (j==i1 || j==i2 || j==i3) continue;	/*--skip on self and directly bonded sites--*/
                	rx = fabs(xi-x[j]);
                	ry = fabs(yi-y[j]);
                	rz = fabs(zi-z[j]);
                	while (rx>Lx/2.0) rx -= Lx;	/*--pbc adjustments--*/
                	while (ry>Ly/2.0) ry -= Ly;
                	while (rz>Lz/2.0) rz -= Lz;
                	r2 = rx*rx + ry*ry + rz*rz;
			if (r2 <= rlist2) { nlist[ptr_val]=j; ptr_val++; }	/*--add site j to list--*/
			if (r2 < 1.0) { /*--If overlap found here: print warning, but keep running--*/
				printf("OVERLAP in list_update()!! [sites %d and %d]\n", i1,j);
				fflush(stdout);
			}
		}
	}
	lst_ptr[N]=ptr_val;	/*--record end of list location--*/
}

void dump_config(ncycle)	/*--Output current system configuration--*/
int ncycle;
{
	int i;

	printf("\nSystem Config after %d MC cycles:\n", ncycle);
	for (i=0; i<NUM; i++)
		printf("%d %f %f %f\n", i, x[i], y[i], z[i]); 
	printf("\n"); fflush(stdout);
}


void gr_setup(rho) /*-- Some initializations for g(r) construction --*/
double rho;	/*--Have removed the (nmol-1)/nmol correction factor to n_unc[]--*/
{
	int i;
	double len, nmol, ai, rgmax, rmin, rmax, Vshell;

	len = (double)LEN; nmol = (double)NMOL;
        for (i=0; i<NG; i++) {	/*--normalization for each g(r) shell--*/
                ai = (double)i;
                rmin = 1.0+ai*DR;
                rmax = rmin+DR;
		Vshell = 4.0*pi*(rmax*rmax*rmax-rmin*rmin*rmin)/3.0;
                /*n_unc[i] = Vshell*len*rho*(nmol-1.0)/nmol;*/
                n_unc[i] = Vshell*len*rho;	/*--uncorrelated fluid result (rho=NMOL/V)--*/
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
        int i1, chn_num, ref_site, j, rindx;
        double xi,yi,zi,rx,ry,rz,r,r2;

        for (i1=0; i1<N; i1++) {	      /*--sweep through all N particles ...--*/
        	xi=x[i1]; yi=y[i1]; zi=z[i1]; /*--(current reference particle)--*/
		chn_num = i1/LEN;		/*--current chain (note use of integer truncation)--*/
		ref_site = chn_num*LEN;		/*--index of 1st site on this chain--*/

        	for (j=0; j<N; j++) {		/*--... and accumulate sums--*/
                	if (j>=ref_site && j<ref_site+LEN) continue; /*--skip on all sites in this molecule--*/
                	rx = fabs(xi-x[j]);
                	ry = fabs(yi-y[j]);
                	rz = fabs(zi-z[j]);
                	while (rx>Lx/2.0) rx -= Lx;	/*--pbc adjustments--*/
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
 
void gr_output(vf)	/*--Construct and output final results--*/
double vf;
{
	int i, nblk;
	double r, diff, g[NBLK][NG], g_avg[NG], g_var[NG];

	/*--zero values--*/
        for (i=0; i<NG; i++) {
		 g_avg[i]=0.0; g_var[i]=0.0;
	}

	/*--Construct averages for each block--*/
	for (nblk=0; nblk<NBLK; nblk++) {
        	for (i=0; i<NG; i++) {	/*--construct g(r)=navg[i]/n_unc[i]--*/
                	g[nblk][i] = n_sum[nblk][i]/(n_unc[i]*(double)(N*ncnt[nblk]));
			g_avg[i] += g[nblk][i];
		}
	}
        for (i=0; i<NG; i++) g_avg[i] /= NBLK;

	/*--Construct variances--*/
        for (i=0; i<NG; i++) g_var[i]=0.0;
	for (nblk=0; nblk<NBLK; nblk++) {
        	for (i=0; i<NG; i++) {
			diff = g_avg[i]-g[nblk][i];
			g_var[i] += diff*diff;
		}
	}
        for (i=0; i<NG; i++) g_var[i]/=(NBLK-1);


        for (i=0; i<NG; i++) {	/*--r-value for g(r)=center of g[i] bin--*/
		r=1.0+((double)i+0.5)*DR;
		if (r*r > rgmax2) break;
                printf("%5.4f  %7.5f +/- %7.5f\n", r, g_avg[i], sqrt(g_var[i]));
	}
}

void chk_bnds()	/*--checks bond lengths for LEN-mer chain fluid--*/
{
	int i, j;
	double dx, dy, dz, r2;

	for (i=0; i<N; i+=LEN) {
		for (j=0; j<NBND; j++) {
			dx=fabs(x[i+j+1]-x[i+j]); dy=fabs(y[i+j+1]-y[i+j]); dz=fabs(z[i+j+1]-z[i+j]);
                	while (dx>Lx/2.0) dx -= Lx;	/*--pbc adjustments--*/
                	while (dy>Ly/2.0) dy -= Ly;
                	while (dz>Lz/2.0) dz -= Lz;
			r2 = dx*dx + dy*dy + dz*dz;
			if (fabs(r2-L2)>1.0e-10) /*--Print warning on detection of incorrect length--*/
				printf("Improper bond length: %d-%d, L=%12.10f\n", i+j, i+j+1, sqrt(r2));
			fflush(stdout);
		}
	}

}

void chkchn_bnds()	/*--checks bond lengths for LEN-mer chain fluid--*/
{
	int i;
	double dx, dy, dz, r2;
	for (i=0; i<NCHN-1; i++) {
        	dx=fabs(x[i+N]-x[i+N+1]); dy=fabs(y[i+N]-y[i+N+1]); dz=fabs(z[i+N]-z[i+N+1]);
               	while (dx>Lx/2.0) dx -= Lx;	/*--pbc adjustments--*/
               	while (dy>Ly/2.0) dy -= Ly;
               	while (dz>Lz/2.0) dz -= Lz;
		r2 = dx*dx + dy*dy + dz*dz;
		if (fabs(r2-L2)>1.0e-10){ /*--Print warning on detection of incorrect length--*/
			printf("Improper bond length: %d-%d, L=%12.10f\n", i+N, i+1+N, sqrt(r2));
			fflush(stdout);
		}
	}
}
