#include "trace_irc.h"

// -------------------------------------------------------------
             /* print help */
void print_help(void){

  printf("\n\n======================================================================\n\n");
  printf("This program traces the intrinsic reaction coordinate from\n");
  printf("the NVAR-dimensional free energy surface generate with a\n");
  printf("CPMD or PLUMED metadynamics simulation, or a 2D potential on a grid\n\n");
  printf("  Requires:\n");
  printf("   1)  colvar_val and dcolvar_val files generaterd with CPMD or\n");
  printf("       HILLS file generaterd with PLUMED or potential on a grid\n");
  printf("   2)  trace_irc.in input file\n");
  printf("  Generates:\n");
  printf("   1)  trace_coursegrid.out, contains the two minima and some points between\n");
  printf("   2)  umbrella.pot, contains all points along the IRC, the free\n");
  printf("                     energy and the forces to use in an CPMD \n");
  printf("                     umbrella sampling run\n");
  printf("   See also the headerfile trace_irc.h for more files and print options\n");
  printf("");
  printf("");
  printf("\n  Keywords required in the inputfile trace_irc:\n");
  printf("    INPUT <TYPE>: Read filenames with potential data from next line(s)\n");
  printf("            TYPE is optional and can be CPMD, PLUMED, or POTENTIAL (default=CPMD)\n");
  printf("            default=colvar_val and dcolvar_val. Old and New CPMD files\n");
  printf("            are recognized. For PLUMED the HILLS file is expected.\n");
  printf("            For POTENTIAL a 2D potential on a regular grid is expected.\n");
  printf("    NVAR:   number of collective variables is read from the next line\n");
  printf("    NGRID:  two integers are read from the next line for the number of\n");
  printf("            grid points. Required if INPUT TYPE is a 2D potential (=POT)\n");
  printf("    MINIMA: next two lines should contain an initial guess for the\n");
  printf("            two minima");
  printf("\n  Additional Keywords recognized in the inputfile trace_irc:\n");
  printf("    SPHERICAL: Spherical hills are being used                 \n");
  printf("    GAUSSCUT: the gaussian cut-off parameter is read from the next\n");
  printf("            line (default: 10^8)\n");
  printf("    BRACKETSIZE: Respectively the bracketsize and the stepsize are\n");
  printf("            read from the next line (used to find the course grid\n");
  printf("            points) (default values: 0.2 0.2). The first number\n");
  printf("            sets a distance right and left from an initial guess\n");
  printf("            between which the energy is minimized for localizing\n");
  printf("            the minima and the irc points. The second number sets\n");
  printf("            the distance between the initial course irc points.\n");
  printf("    NSMALLSTEPS: Number of points to generate in between the course\n");
  printf("            grid that is traced first (default 10 points). The minimum\n");
  printf("            is 1 small step per segment, taking just the coarse grid points.\n");
  printf("    NSIDEPOINTS: Number of points to generate outside the two minima\n");
  printf("            with the fine grid spacing (default 10 points)\n");
  printf("            (NB this spacing equals the stepsize (see BRACKETSIZE) or the\n");
  printf("            distance between the minima (which ever smaller) divided by \n");
  printf("            NSMALLSTEPS, but the final distance between the found points\n");
  printf("            is NOT fixed (only their projection on the CG vectors is)).\n");
  printf("    TEDGECONST: if set then the NSIDEPOINTS have the same energy as the\n");
  printf("                minima so that the force will vanish outside the minima\n");
  printf("    SMOOTHFACTOR: 2 gaussian widths are read for smoothing respec-\n");
  printf("            tively the points along the IRC and the energy (required\n");
  printf("            to compute smooth forces (default no smoothing)\n");
  printf("    TNOMINOPT: if set then no optimization of the minima \n");
  printf("    TNOTRACE: if set then locate minima only\n");
  printf("    TMONTECARLO: performs the integration of the free energy surface\n");
  printf("             over the nvar-1 directions perpendicular to the IRC\n");
  printf("             and an effective umbrella potential is computed instead\n");
  printf("             of just the potential at the minimum path (IRC). Reads 1\n");
  printf("             integer and 1 float from the next line for the number of \n");
  printf("             accepted MC moves and the maximum stepsize\n");
  printf("     FES_MOVIE: Interval in number of hills to write the evolution of the\n");
  printf("             free energy profile along the final path\n");
  printf("");
  printf("");
  printf("");
  printf("\n\n======================================================================\n\n");

}

// -------------------------------------------------------------

double myrandom(idumm){

  double random;
  const double d2p31m=2147483647.0, d2p31=2147483648.0, d7p5=16807.0;

  dseed=fmod(d7p5*dseed,d2p31m);
  random=dseed/d2p31;
  return random;

}

// -------------------------------------------------------------
                      /* radius */
double radius(double *x, double *y, int ndim){
  int i;
  double r=0.0, dx;

  for(i=0;i<ndim;i++){
    dx=x[i]-y[i];
    r=r+dx*dx;
  }
  return sqrt(r);
}

// -------------------------------------------------------------
                      /* read input */
int read_input(HILLSPATH *hills, double **points, double *bracketsize, double *stepsize, int *nsteps, double *gausswidth, int *nside, int *tedgeconst, int *tminopt, int *ttrace, int *tmc, int *maxmcmoves, double *maxstep, double *drmax, double *mcminetot, char *colvarfile, char *hillfile,int *inputtype,int *nprint_fesmovie){

  FILE *fpin;
  char keyword[50], line[MAXLINE];
  int  i, nvar;

  nvar=hills->nvar;
  if( (fpin = fopen(INPUTFILE,"r")) == NULL){
    printf("Cannot open %s \nStopped\n",INPUTFILE);
    print_help();
    return EXIT_FAILURE;
  }
  while( fgets(line,MAXLINE,fpin) != NULL ){
    sscanf(strtok(line," "),"%s",keyword);
    //    printf("keyword: %s\n",keyword);
    if( strcasecmp(keyword,"NVAR") == 0 ){
      fscanf(fpin,"%i",&nvar);
      printf("  Number of CVs:             %d\n",nvar);
      fgets(line,MAXLINE,fpin);
    }
    else if( strcasecmp(keyword,"MINIMA") == 0 ){
      if(nvar==-1){
        printf(" ERROR: Give NVAR before MINIMA in %s\n Stopped\n",INPUTFILE);
        return EXIT_FAILURE;
      }
      fgets(line,MAXLINE,fpin);
      sscanf(strtok(line," "),"%lf",&points[0][0]);
      for(i=1;i<nvar;i++) sscanf(strtok(NULL," "),"%lf",&points[0][i]);
      printf("  Guess minimim 1:           %lg",points[0][0]);
      for(i=1;i<nvar;i++) printf(", %lg",points[0][i]); printf("\n");
      fgets(line,MAXLINE,fpin);
      sscanf(strtok(line," "),"%lf",&points[1][0]);
      for(i=1;i<nvar;i++) sscanf(strtok(NULL," "),"%lf",&points[1][i]);
      printf("  Guess minimim 2:           %lg",points[1][0]);
      for(i=1;i<nvar;i++) printf(", %lg",points[1][i]); printf("\n");
    }
    else if( strcasecmp(keyword,"INPUT") == 0 ){
      if(sscanf(strtok(NULL," "),"%s",keyword) != 1){
	*inputtype=0;  /* default is CPMD input */
	fscanf(fpin,"%s%s",colvarfile,hillfile);
	printf(" Using CPMD input from %s and %s\n",colvarfile,hillfile);
      }else{
	if(strcasecmp(keyword,"CPMD") == 0){
	  *inputtype=0;
	  fscanf(fpin,"%s%s",colvarfile,hillfile);
	  printf(" Using CPMD input from %s and %s\n",colvarfile,hillfile);
	}else if(strcasecmp(keyword,"PLUMED") == 0){
	  *inputtype=1;
	  fscanf(fpin,"%s",hillfile);
	  printf("  Using PLUMED input from %s\n",hillfile);
	}else if(strcasecmp(keyword,"POTENTIAL") == 0){
	  *inputtype=2;
	  fscanf(fpin,"%s",hillfile);
	  printf("  Using potential input from %s\n",hillfile);
	}else{
	  printf(" ERROR: INPUT %s unknown (use CPMD, PLUMED, or POTENTIAL)\n",keyword);
	  return EXIT_FAILURE;
	}
      }
      fgets(line,MAXLINE,fpin);
    }
    else if( strcasecmp(keyword,"GAUSSCUT") == 0 ){
      fscanf(fpin,"%lf",&hills->cutgauss);
      fgets(line,MAXLINE,fpin);
    }
    else if( strcasecmp(keyword,"SPHERICAL") == 0 ){
      hills->tspherical=1;
    }
    else if( strcasecmp(keyword,"SYMMETRY") == 0 ){
      hills->tsym=1;
      if( sscanf(&line[9],"%d",&hills->isym) != 1 ){
	hills->isym=-1;
	printf("  All symmetry equivalent hills are used\n");
      }else{
	printf("  Only hills of symmetry point %d are used\n",hills->isym);
      }
    }
    else if( strcasecmp(keyword,"BRACKETSIZE") == 0 ){
      fscanf(fpin,"%lf%lf",bracketsize,stepsize);
      printf("  Bracket and step size:     %lf %lf\n",*bracketsize,*stepsize);
      fgets(line,MAXLINE,fpin);
    }
    else if( strcasecmp(keyword,"NSMALLSTEPS") == 0 ){
      fscanf(fpin,"%d",nsteps);
      if( *nsteps < 1) *nsteps=1;
      printf("  Number of small steps:     %d\n",*nsteps);
      fgets(line,MAXLINE,fpin);
    }
    else if( strcasecmp(keyword,"NSIDEPOINTS") == 0 ){
      fscanf(fpin,"%d",nside);
      printf("  Number of flanking points: %d\n",*nside);
      fgets(line,MAXLINE,fpin);
    }
    else if( strcasecmp(keyword,"TEDGECONST") == 0 ){
      *tedgeconst=1;
    }
    else if( strcasecmp(keyword,"SMOOTHFACTOR") == 0 ){
      fscanf(fpin,"%lf%lf",&gausswidth[0],&gausswidth[1]);
      printf("  Smooth factors:            %lf %lf\n",gausswidth[0],gausswidth[1]);
      fgets(line,MAXLINE,fpin);
    }
    else if( strcasecmp(keyword,"TNOMINOPT") == 0 ){
      *tminopt=0;
    }
    else if( strcasecmp(keyword,"TNOTRACE") == 0 ){
      *ttrace=0;
    }
    else if( strcasecmp(keyword,"TMONTECARLO") == 0 ){
      *tmc=1;
      fscanf(fpin,"%d%lf",maxmcmoves,maxstep);
      fgets(line,MAXLINE,fpin);
    }
    else if( strcasecmp(keyword,"MCMAXVALLEYWIDTH") == 0 ){
      fscanf(fpin,"%lf",drmax);
      fgets(line,MAXLINE,fpin);
    }
    else if( strcasecmp(keyword,"MCMINETOT") == 0 ){
      fscanf(fpin,"%lf",mcminetot);
      fgets(line,MAXLINE,fpin);
    }
    else if( strcasecmp(keyword,"NEWCODE") == 0 ){
      hills->tnewcode=1;
    }
    else if( strcasecmp(keyword,"FES_MOVIE") == 0 ){
      fscanf(fpin,"%d",nprint_fesmovie);
      fgets(line,MAXLINE,fpin);
    }
    else if( strcasecmp(keyword,"NGRID") == 0 ){
      fgets(line,MAXLINE,fpin);
      if( sscanf(line,"%d%d",&hills->nx1,&hills->nx2) != 2){
	printf("\n  ERROR: failed to read two integers after the NGRID keyword.\n\n");
	exit(1);
      }
      printf("  Potential grid size:       %d x %d\n",hills->nx1,hills->nx2);
    }
    else{
      printf("  Non-fatal warning: Keyword not recognized: %s \n",keyword);
	//	return EXIT_FAILURE;
    }
  } 
  hills->nvar = nvar;
  fclose(fpin);

  /* Sanity checks */
  if( (*inputtype==2) && (hills->nx2==-1) ){
    printf("\n  ERROR: input type POTENTIAL selected but no NGRID specified.\n\n");
    exit(1);
  }


  /* Some more feedback */
  if(hills->tspherical){
    printf("  Spherical hills:           yes\n");
  }else{
    printf("  Spherical hills:           no\n");
  }
  if(*tedgeconst){
    printf("  Sides constant energy:     yes\n");
  }else{
    printf("  Sides constant energy:     no\n");
  }
  if(gausswidth[0]==-999.){
    printf("  Smooting of path:          no\n");
  }
  if(gausswidth[1]==-999.){
    printf("  Smooting of energy:        no\n");
  }
  if(*tminopt){
    printf("  Optimize minima:           yes\n");
  }else{
    printf("  Optimize minima:           no\n");
  }
  if(*ttrace){
    printf("  Trace path:                yes\n");
  }else{
    printf("  Trace path:                no\n");
  }
  if(*nprint_fesmovie > 0){
    printf("  Write FES with interval of %d HILLS\n",*nprint_fesmovie);
  }
  if(*tmc){
    printf("  MC integration:            yes\n");
    printf("  Max moves and stepsize:    %d, %lf",*maxmcmoves,*maxstep);
  }else{
    printf("  MC integration:            no\n");
  }
  return EXIT_SUCCESS;
}

// -------------------------------------------------------------
  /* open colvar_val and dcolvar_val to read the hills */
int read_hills_from_CPMDfile(HILLSPATH *hills, char *colvarfile,char *hillfile){

  FILE *fpcolvar, *fpdcolvar;
  char line[MAXLINE];
  int i, ihill, istep, nvar, ndummy,isym;
  double dummy;

  nvar=hills->nvar;
  if( (fpcolvar = fopen(colvarfile,"r")) == NULL){
    printf("\n  ERROR: Failed to open colvar file %s \nStopped\n\n",colvarfile);
    return EXIT_FAILURE;
  }
  if( (fpdcolvar = fopen(hillfile,"r")) == NULL){
    printf("\n  ERROR: Failed to open hills file %s \nStopped\n\n",hillfile);
    return EXIT_FAILURE;
  }

  /* allocate memory for the colvar trajectory and scale factors and read them in */
  ihill=0;
  hills->cvpath = (double **) malloc(nvar*sizeof(double *));
  for(i=0;i<nvar;i++) hills->cvpath[i] = (double *) malloc(NHILLSMAX*sizeof(double));
  hills->scale = (double **) malloc(nvar*sizeof(double *));
  for(i=0;i<nvar;i++) hills->scale[i] = (double *) malloc(NHILLSMAX*sizeof(double));
  while( fgets(line,MAXLINE,fpcolvar) != NULL ){
    sscanf(strtok(line," "),"%d",&istep);
    if(hills->tsym) sscanf(strtok(NULL," "),"%d",&isym);
    if( (hills->isym==-1) || (hills->isym==isym) ){
      for(i=0;i<nvar;i++){
        sscanf(strtok(NULL," "),"%lg",&hills->cvpath[i][ihill]);
      }
      for(i=0;i<nvar;i++){
        sscanf(strtok(NULL," "),"%lg",&hills->scale[i][ihill]);
      }
      ihill++;
      if(ihill >= NHILLSMAX){
        printf("Number of hills to read in more than %d\n",NHILLSMAX);
        printf("Increase MAXNHILLS\nProgram stopped\n");
        return EXIT_FAILURE;
      }
    }
  }  
  hills->nhills=ihill;
  printf("  Number of hills:           %d\n",hills->nhills);

  /* allocate memory for the hill height and width and read them in */
  if(hills->tnewcode){
    ndummy=1;
  } else{
    ndummy=nvar+1;
  }
  hills->hllw = (double *) malloc(hills->nhills*sizeof(double));
  hills->hllh = (double *) malloc(hills->nhills*sizeof(double));
  for(ihill=0;ihill<hills->nhills;){
    fgets(line,MAXLINE,fpdcolvar);
    sscanf(strtok(line," "),"%d",&istep);
    if(hills->tsym) sscanf(strtok(NULL," "),"%d",&isym);
    if( (hills->isym==-1) || (hills->isym==isym) ){    
      for(i=0;i<ndummy;++i) sscanf(strtok(NULL," "),"%lg",&dummy);
      sscanf(strtok(NULL," "),"%lg",&hills->hllw[ihill]);
      sscanf(strtok(NULL," "),"%lg",&hills->hllh[ihill]);
      ihill++;
    }
  }
  fclose(fpdcolvar);
  fclose(fpcolvar);

  printf("  Hill width and heigth of last hill are %lf, %lf\n",hills->hllw[ihill-1],hills->hllh[ihill-1]);
  
  return EXIT_SUCCESS;
}

// -------------------------------------------------------------
  /* open PLUMED HILLS file to read the hills */
int read_hills_from_PLUMEDfile(HILLSPATH *hills,char *hillfile){

  FILE *fphills;
  char line[MAXLINE], *pt;
  int i, ihill, nvar, ndummy,isym, nhills;
  double dummy, time;

  nvar=hills->nvar;
  if( (fphills = fopen(hillfile,"r")) == NULL){
    printf("\n  ERROR: Failed to open hills file %s \nStopped\n\n",hillfile);
    return EXIT_FAILURE;
  }

  /* count number of hills */
  nhills=0;
  while( fgets(line,MAXLINE,fphills) != NULL ){
    pt=strtok(line," ");   /* skip comment lines starting with # */
    if( pt[0] !='#' ) nhills++;
  }
  //  printf("  Number of hills:           %d\n",nhills);
  rewind(fphills);

  /* allocate memory for the hills trajectory and scale factors */
  ihill=0;
  hills->cvpath = (double **) malloc(nvar*sizeof(double *));
  for(i=0;i<nvar;i++) hills->cvpath[i] = (double *) malloc(nhills*sizeof(double));
  hills->scale = (double **) malloc(nvar*sizeof(double *));
  for(i=0;i<nvar;i++) hills->scale[i] = (double *) malloc(nhills*sizeof(double));
  hills->hllw = (double *) malloc(nhills*sizeof(double));
  hills->hllh = (double *) malloc(nhills*sizeof(double));

  /* read hills */
  while( fgets(line,MAXLINE,fphills) != NULL ){
    pt=strtok(line," ");
    if( pt[0] !='#' ){   /* skip comment lines starting with # */
      sscanf(pt,"%lg",&time);
      if(hills->tsym) sscanf(strtok(NULL," "),"%d",&isym);
      if( (hills->isym==-1) || (hills->isym==isym) ){
	for(i=0;i<nvar;i++){
	  sscanf(strtok(NULL," "),"%lg",&hills->cvpath[i][ihill]);
	}
	for(i=0;i<nvar;i++){
	  sscanf(strtok(NULL," "),"%lg",&hills->scale[i][ihill]);
	}
	sscanf(strtok(NULL," "),"%lg",&hills->hllh[ihill]);
	hills->hllw[ihill] = 1.0;
	ihill++;
      }
    }
  }  
  hills->nhills=ihill;
  printf("  Number of hills:           %d\n",hills->nhills);

  fclose(fphills);

  printf("  Hill heigth and widths of last hill are %lf and %lf",hills->hllh[ihill-1],hills->scale[0][ihill-1]);
  for(i=1;i<nvar;i++) printf(", %lf",hills->scale[i][ihill-1]);
  printf("\n");
  
  return EXIT_SUCCESS;
}

// -------------------------------------------------------------
/* open POTENTIAL file and read input data */
int read_2Dpotential(HILLSPATH *hills,char *potfile){

  int i, ix, iy, nx1, nx2;
  double dx1, dx2;
  FILE *fpin;

  /* allocate memory */
  nx1=hills->nx1;
  nx2=hills->nx2;
  hills->x1pot = (double *)malloc(nx1*sizeof(double));
  hills->x2pot = (double *)malloc(nx2*sizeof(double));
  hills->extpot = (double **)malloc(nx1*sizeof(double *));
  hills->extdvxx = (double **)malloc(nx1*sizeof(double *));
  hills->extdvx1 = (double **)malloc(nx1*sizeof(double *));
  hills->extdvx2 = (double **)malloc(nx1*sizeof(double *));
  hills->extpot[0] = (double *)malloc(nx1*nx2*sizeof(double));
  hills->extdvxx[0] = (double *)malloc(nx1*nx2*sizeof(double));
  hills->extdvx1[0] = (double *)malloc(nx1*nx2*sizeof(double));
  hills->extdvx2[0] = (double *)malloc(nx1*nx2*sizeof(double));
  for(i=1;i<nx1;i++){
    hills->extpot[i] = hills->extpot[i-1] + nx2;
    hills->extdvxx[i] = hills->extdvxx[i-1] + nx2;
    hills->extdvx1[i] = hills->extdvx1[i-1] + nx2;
    hills->extdvx2[i] = hills->extdvx2[i-1] + nx2;
  }

  /* open file */
  if( (fpin = fopen(potfile,"r")) == NULL){
    printf("\n  ERROR: Failed to open potential file %s \nStopped\n\n",potfile);
    return EXIT_FAILURE;
  } 

  /* read first line and get grid origin */
  if( fscanf(fpin,"%lf%lf%lf",&hills->x1pot[0],&hills->x2pot[0],&hills->extpot[0][0]) != 3){
    printf("\n  ERROR: Failed to read first line from potential file %s\n\n",potfile);
    exit(1);
  }
  rewind(fpin);

  /* read file */
  for(ix=0;ix<nx1;ix++){
    for(iy=0;iy<nx2;iy++){
      if( fscanf(fpin,"%lf%lf%lf",&hills->x1pot[ix],&hills->x2pot[iy],&hills->extpot[ix][iy]) != 3){
	printf("\n  ERROR: Failed while reading line %d from potential file %s\n\n",ix*iy,potfile);
	exit(1);
      }
    }
  }

  /* negate the potential (contrary to metadynamics Gaussians, external potentials are typically negative valued) */
  for(ix=0;ix<nx1;ix++){
    for(iy=0;iy<nx2;iy++){
      hills->extpot[ix][iy] *= -1.0;
    }
  }  

  /* some feedback */
  dx1 = (hills->x1pot[nx1-1]-hills->x1pot[0])/(nx1-1);
  dx2 = (hills->x2pot[nx2-1]-hills->x2pot[0])/(nx2-1);
  printf("  Potential grid origin:     %lf %lf\n",hills->x1pot[0],hills->x2pot[0]);
  printf("  Potential grid interval:   %lf %lf\n",dx1,dx2);

  fclose(fpin);
  return EXIT_SUCCESS;
}

// -------------------------------------------------------------
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
double **dmatrix(int nrl,int nrh,int ncl,int nch)
{
  int i,nrow,ncol;
  double **m;

  nrow = nrh-nrl+1;  ncol = nch-ncl+1;

  if(nrow*ncol==0) return (double **)NULL;
  /* allocate pointers to rows */

  m = (double **)malloc(nrow*sizeof(double*));
  if(!m){
    fprintf(stderr,"ERROR: can't allocate memory for matrix (%lu bytes)\n",
            nrow*sizeof(double*));
    exit(1);
  }
  m -= nrl;

  /* allocate rows */
  m[nrl] = (double *)calloc(nrow*ncol,sizeof(double));
  if(!m[nrl]){
    fprintf(stderr,
            "ERROR: can't allocate memory for arrays for matrix (%lu bytes)\n",
            nrow*ncol*sizeof(double));
    exit(1);
  }
  m[nrl] -= ncl;

  /* set pointers to rows */
  for(i=nrl+1;i<=nrh;i++) m[i] = m[i-1]+ncol;

  return m;
}

// -------------------------------------------------------------
void free_dmatrix(double **m,int nrl,int nrh,int ncl,int nch)
{
  free((char *)(m[nrl]+ncl)); free((char *)(m+nrl));
}

// -------------------------------------------------------------
void bcucof(double y[], double y1[], double y2[], double y12[], double d1, double d2,
	    double **c){

  static int wt[16][16]=
    { 1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
      -3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0,
      2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0,
      0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
      0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1,
      0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1,
      -3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0,
      9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2,
      -6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2,
      2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0,
      -6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1,
      4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1};
  int l,k,j,i;
  double xx,d1d2,cl[16],x[16];

  d1d2=d1*d2;
  for (i=1;i<=4;i++) {
    x[i-1]=y[i];
    x[i+3]=y1[i]*d1;
    x[i+7]=y2[i]*d2;
    x[i+11]=y12[i]*d1d2;
  }
  for (i=0;i<=15;i++) {
    xx=0.0;
    for (k=0;k<=15;k++) xx += wt[i][k]*x[k];
    cl[i]=xx;
  }
  l=0;
  for (i=1;i<=4;i++)
    for (j=1;j<=4;j++) c[i][j]=cl[l++];
}

// -------------------------------------------------------------
void bcuint(double y[], double y1[], double y2[], double y12[], double x1l,
	    double x1u, double x2l, double x2u, double x1, double x2, double *ansy,
	    double *ansy1, double *ansy2){

  void bcucof(double y[], double y1[], double y2[], double y12[], double d1,
	      double d2, double **c);
  int i;
  double t,u,d1,d2,**c;
  
  c=dmatrix(1,4,1,4);
  d1=x1u-x1l;
  d2=x2u-x2l;
  bcucof(y,y1,y2,y12,d1,d2,c);
  if (x1u == x1l || x2u == x2l) {
    fprintf(stderr,"Bad input in routine bcuint");
    exit(1);
  } 
  t=(x1-x1l)/d1;
  u=(x2-x2l)/d2;
  *ansy=(*ansy2)=(*ansy1)=0.0;
  for (i=4;i>=1;i--) {
    *ansy=t*(*ansy)+((c[i][4]*u+c[i][3])*u+c[i][2])*u+c[i][1];
    *ansy2=t*(*ansy2)+(3.0*c[i][4]*u+2.0*c[i][3])*u+c[i][2];
    *ansy1=u*(*ansy1)+(3.0*c[4][i]*t+2.0*c[3][i])*t+c[2][i];
  }
  *ansy1 /= d1;
  *ansy2 /= d2;
  free_dmatrix(c,1,4,1,4);
}

/*--------------------------------------------------------------*/
void get_2Dextpot_derivatives(HILLSPATH *hills){

  int i, j;
  double **extdvx1 = hills->extdvx1;
  double **extdvx2 = hills->extdvx2;
  double **extdvxx = hills->extdvxx;
  double **extpot = hills->extpot;
  double *x1pot = hills->x1pot;
  double *x2pot = hills->x2pot;

  for(i=0;i<hills->nx1;i++){
    extdvx1[i][0] =  extdvx1[i][hills->nx2-1] = 0.;
    extdvx2[i][0] =  extdvx2[i][hills->nx2-1] = 0.;
    extdvxx[i][0] =  extdvxx[i][hills->nx2-1] = 0.;
  }
  for(j=0;j<hills->nx2;j++){
    extdvx1[0][j] =  extdvx1[hills->nx1-1][j] = 0.;
    extdvx2[0][j] =  extdvx2[hills->nx1-1][j] = 0.;
    extdvxx[0][j] =  extdvxx[hills->nx1-1][j] = 0.;
  }
  for(i=1;i<hills->nx1-1;i++){
    for(j=1;j<hills->nx2-1;j++){
      extdvx1[i][j] = (extpot[i+1][j] - extpot[i-1][j]) / (x1pot[i+1] - x1pot[i-1]);
      extdvx2[i][j] = (extpot[i][j+1] - extpot[i][j-1]) / (x2pot[j+1] - x2pot[j-1]);
      extdvxx[i][j] = (extpot[i+1][j+1] - extpot[i+1][j-1] - extpot[i-1][j+1] + extpot[i-1][j-1])
        / ((x1pot[i+1] - x1pot[i-1]) * (x2pot[j+1] - x2pot[j-1]));
    }
  }
}

// -------------------------------------------------------------
double energy_pot(double *point, HILLSPATH *hills){

  int i, j;
  double y[5], y1[5], y2[5], y3[5], v, vdx1, vdx2;
  double *x1pot = hills->x1pot;
  double *x2pot = hills->x2pot;
  double **extpot = hills->extpot;
  double **extdvx1 = hills->extdvx1;
  double **extdvx2 = hills->extdvx2;
  double **extdvxx = hills->extdvxx;

  /* find in which grid square the point is */
  i = j = 0;
  if((point[0] < x1pot[0]) || (point[0] > x1pot[hills->nx1-1]) || (point[1] < x2pot[0]) || (point[1] > x2pot[hills->nx2-1])){
    /* pcolvar is off the potential grid where the external force is zero */
    return 0.0;
  }

  i = j = 1;
  while(point[0] > x1pot[i]) i++;
  while(point[1] > x2pot[j]) j++;
  i--;
  j--;
     
  y[1] = extpot[i][j];
  y[2] = extpot[i+1][j];
  y[3] = extpot[i+1][j+1];
  y[4] = extpot[i][j+1];

  y1[1] = extdvx1[i][j];
  y1[2] = extdvx1[i+1][j];
  y1[3] = extdvx1[i+1][j+1];
  y1[4] = extdvx1[i][j+1];

  y2[1] = extdvx2[i][j];
  y2[2] = extdvx2[i+1][j];
  y2[3] = extdvx2[i+1][j+1];
  y2[4] = extdvx2[i][j+1];

  y3[1] = extdvxx[i][j];
  y3[2] = extdvxx[i+1][j];
  y3[3] = extdvxx[i+1][j+1];
  y3[4] = extdvxx[i][j+1];

  bcuint(y,y1,y2,y3,x1pot[i],x1pot[i+1],x2pot[j],x2pot[j+1],
    point[0],point[1],&v,&vdx1,&vdx2);

#ifdef DEBUG
  printf("punt %d (%d,%d): x=%lf y=%lf E=%lf (%lf %lf %lf)\n",1,i,j,x1pot[i],x2pot[j],y[1],y1[1],y2[1],y3[1]);
  printf("punt %d (%d,%d): x=%lf y=%lf E=%lf (%lf %lf %lf)\n",2,i+1,j,x1pot[i+1],x2pot[j],y[2],y1[2],y2[2],y3[2]);
  printf("punt %d (%d,%d): x=%lf y=%lf E=%lf (%lf %lf %lf)\n",3,i+1,j+1,x1pot[i+1],x2pot[j+1],y[3],y1[3],y2[3],y3[3]);
  printf("punt %d (%d,%d): x=%lf y=%lf E=%lf (%lf %lf %lf)\n",4,i,j+1,x1pot[i],x2pot[j+1],y[4],y1[4],y2[4],y3[4]);
  printf("Potential at (%lf,%lf) = %g; Forces are %g, %g\n",point[0],point[1],v,vdx1,vdx2);
#endif

  //  colvar->fcolvar[0] -= vdx1;
  //  colvar->fcolvar[1] -= vdx2;

  return v;

}

// -------------------------------------------------------------
double energy_meta(double *coord, HILLSPATH *hills){

  double gperp, gperp_cut, gpar, expn, expn2, sigma2;
  double diffpos[NVARMAX], deltacv[NVARMAX], stepgauss;
  int i, ihill, nhills;

  nhills=hills->nhills;
  stepgauss=0.0;
  if(hills->tspherical==1){
    for(ihill=0;ihill<nhills-1;ihill++){
      expn=0.0;
      for(i=0;i<hills->nvar;i++){
	diffpos[i]=(coord[i]-hills->cvpath[i][ihill])/hills->scale[i][ihill];
	expn=expn+diffpos[i]*diffpos[i];
      }
      expn=sqrt(expn);
      if(expn < hills->hllw[i]*hills->cutgauss){
	gperp=exp(-0.5*pow(expn/hills->hllw[ihill],2.0));
	gperp_cut=exp(-0.5*pow(hills->cutgauss,2.0));
	stepgauss += hills->hllh[ihill]*(gperp-gperp_cut);
      }
    }
  } else{
    for(ihill=0;ihill<nhills-1;ihill++){
      expn=0.0;
      expn2=0.0;
      sigma2=0.0;
      for(i=0;i<hills->nvar;i++){
	diffpos[i]=(coord[i]-hills->cvpath[i][ihill])/hills->scale[i][ihill];
	deltacv[i]=(hills->cvpath[i][ihill+1]-hills->cvpath[i][ihill])/hills->scale[i][ihill];
	expn=expn+diffpos[i]*diffpos[i];
	expn2=expn2+diffpos[i]*deltacv[i];
	sigma2=sigma2+deltacv[i]*deltacv[i];
      }
      expn=sqrt(expn);
      sigma2=sigma2*sigma2;
      if(sigma2 < 1.0e-6) sigma2=1.0e-6;
      if(expn < hills->hllw[i]*hills->cutgauss){
	gperp=exp(-0.5*pow(expn/hills->hllw[ihill],2.0));
	gperp_cut=exp(-0.5*pow(hills->cutgauss,2.0));
	gpar=hills->hllh[ihill]*exp(-0.5*(expn2*expn2/sigma2));
	stepgauss += gpar*(gperp-gperp_cut);
      }
    }
  }
  return stepgauss;
}

// -------------------------------------------------------------
void print_point(FILE *fpout, HILLSPATH *hills, double *p){
  int i;

  fprintf(fpout,PRINTFORMAT,p[0]);
  for(i=1;i<hills->nvar;i++) fprintf(fpout,PRINTFORMAT,p[i]);
  fprintf(fpout,PRINTFORMAT,energy(p,hills));
  fprintf(fpout,"\n");
  fflush(fpout);
}

// -------------------------------------------------------------
void print_movie(HILLSPATH *hills, double **irc_traj,int totpoints, int nprint){
  int i, ivar, nvar, ipoint, nhillstot;
  char filename[MAXLINE];
  FILE *fpout;

  nhillstot = hills->nhills;    /* temporary store the actual total number of hills */
  hills->nhills = nprint;
  nvar = hills->nvar;

  while(hills->nhills <= nhillstot){
    sprintf(filename,"fes_%d.out",hills->nhills);
    printf("  Writing energy profile along path to file %s. Number of hills: %d\n",filename,hills->nhills);

    fpout = fopen(filename,"w");

    for(ipoint=1;ipoint<totpoints-1;ipoint++){        /* force in points between minima */
      fprintf(fpout,"  %d",ipoint+1);
      for(ivar=0;ivar<nvar;ivar++) fprintf(fpout,"   %lf",irc_traj[ipoint][ivar]);
      fprintf(fpout,"   %lf\n",energy(irc_traj[ipoint],hills));
    }
    fclose(fpout);
    hills->nhills += nprint;
  }    

  hills->nhills = nhillstot;   /* reset to the actual total number of hills */
  return;
}
 
// -------------------------------------------------------------
int bracket_minimum(HILLSPATH *hills, double *min, double bracketsize){

  double pointright[NVARMAX], pointleft[NVARMAX];
  double Eright, Eleft, Emiddle;
  int  i,  ivar, nvar, iconsistency;

  nvar=hills->nvar;
  iconsistency=1;
  printf("============= Bracketing minimum ===============================\n");
  printf("ivar  left bracket E       minimum Energy       right bracket E \n");
  while(iconsistency){
    iconsistency=0;
    ivar=0;
    while(ivar<nvar){
      for(i=0;i<nvar;i++)pointleft[i]=min[i] ;
      for(i=0;i<nvar;i++)pointright[i]=min[i] ;
      pointright[ivar] = pointright[ivar] + bracketsize;
      pointleft[ivar] = pointleft[ivar] - bracketsize;
      Eleft=energy(pointleft,hills);
      Eright=energy(pointright,hills);
      Emiddle=energy(min,hills);
      printf("  %d  %lf %lf    %lf %lf    %lf %lf\n",ivar,pointleft[ivar],Eleft,min[ivar],Emiddle,pointright[ivar],Eright);
      if((Eleft>Emiddle)&&(Eright>Emiddle)){
	printf("Bracketing failed.\nAdjust bracketsize.\nProgram stopped.\n");
	return EXIT_FAILURE;
      }
      if(Eleft>Emiddle){
	min[ivar] = pointleft[ivar];
	iconsistency=1;
      }
      else if(Eright>Emiddle){
	min[ivar] = pointright[ivar];
	iconsistency=1;
      }
      else{
	ivar++;
      }
    }
  }
  return EXIT_SUCCESS;
}

// -------------------------------------------------------------
void minimize_energy(HILLSPATH *hills, double *min, double bracketsize){

  double pointright[NVARMAX], pointleft[NVARMAX], point[NVARMAX];
  double step, Eright, Eleft, Emiddle, Eprevious;
  int  i, ihill, ivar, icycle, icycle2, nvar;

  nvar=hills->nvar;
  icycle2=0;
  printf("============= Minimizing =======================================\n");
  printf("cycle   var.1 ");
  for(ivar=1;ivar<nvar;ivar++) printf("   var.%d ",ivar+1);
  printf("   Energy \n");  
  printf("> %d   %lf",icycle2,min[0]);
  for(i=1;i<nvar;i++) printf(" %lf",min[i]);
  printf("  %lg\n",energy(min,hills));
  Emiddle=0.0;
  Eprevious=Emiddle-2.0*ETOL;
  while( ((Emiddle-Eprevious) > ETOL) && (icycle2<100)){
    Eprevious=Emiddle;
    for(ivar=0;ivar<nvar;ivar++){
      for(i=0;i<nvar;i++) pointleft[i] = min[i] ;
      for(i=0;i<nvar;i++) pointright[i] = min[i] ;
      for(i=0;i<nvar;i++) point[i] = min[i] ;
      pointright[ivar] = pointright[ivar] + bracketsize;
      pointleft[ivar] = pointleft[ivar] - bracketsize;
      icycle=0;
      while( (pointright[ivar]-pointleft[ivar]) > TOL){
	Eleft=energy(pointleft,hills);
	Eright=energy(pointright,hills);
	Emiddle=energy(min,hills);
	// 	printf("%d %d  %d  %lf %lf  %lf %lf  %lf %lf\n",icycle2,ivar,icycle,pointleft[ivar],Eleft,min[ivar],Emiddle,pointright[ivar],Eright);
	if(Eleft>Eright){
	  step=0.5*(pointright[ivar] - min[ivar]);
	  point[ivar] = pointright[ivar] - step;
	  if(energy(point,hills)<Emiddle){
	    pointright[ivar]=point[ivar];
	  }
	  else{
	    pointleft[ivar]=min[ivar];
	    min[ivar]=point[ivar];
	  } 
	}
	else{
	  step=0.5*(min[ivar]-pointleft[ivar]);
	  point[ivar] = pointleft[ivar] + step;
	  if(energy(point,hills)<Emiddle){
	    pointleft[ivar]=point[ivar];
	  }
	  else{
	    pointright[ivar]=min[ivar];
	    min[ivar]=point[ivar];
	  } 
	}
	icycle++;
      }
      /*    printf("%d  %lf",icycle,point[0]);
	    for(i=1;i<nvar;i++) printf("  %lf",point[i]);
	    printf("  %lg\n",energy(point,hills));    */
    }
    printf("> %d   %lf",icycle2+1,min[0]);
    for(i=1;i<nvar;i++) printf(" %lf",min[i]);
    printf("  %lg\n",energy(min,hills));
    icycle2++;
  }      
}

// -------------------------------------------------------------
int bracket_minimum_onsphere(HILLSPATH *hills, double *min, double* center, double stepsize, double bracketsize, int ivar ,FILE *fpout){

  double pright[NVARMAX], pleft[NVARMAX];
  double Eright, Eleft, Emiddle, fact;
  int  i, nvar, icycle;

  nvar=hills->nvar;
  icycle=0;
  //  printf("============= Bracketing minimum on sphere =====================\n");
  //  printf("ivar  left bracket E       minimum Energy       right bracket E \n");
  while(icycle<MAXCYCLE){
    icycle++;
    /* make right bracket point */
    for(i=0;i<nvar;i++) pright[i]=min[i] ;
    pright[ivar] = pright[ivar] + bracketsize;
    fact=stepsize/radius(pright,center,nvar);
    for(i=0;i<nvar;i++) pright[i]= center[i]+(pright[i]-center[i])*fact;
    //    print_point(fpout,hills,pright);

    /* make left bracket point */
    for(i=0;i<nvar;i++) pleft[i]=min[i] ;
    pleft[ivar] = pleft[ivar] - bracketsize;
    //    print_point(fpout,hills,pleft);
    fact=stepsize/radius(pleft,center,nvar);
    for(i=0;i<nvar;i++) pleft[i]= center[i]+(pleft[i]-center[i])*fact;
    //    print_point(fpout,hills,pleft);

    Eleft=energy(pleft,hills);
    Eright=energy(pright,hills);
    Emiddle=energy(min,hills);
    //    printf("  %d  %lf %lf    %lf %lf    %lf %lf\n",ivar,pleft[ivar],Eleft,min[ivar],Emiddle,pright[ivar],Eright);
    if((Eleft<Emiddle)&&(Eright<Emiddle)){
	return EXIT_SUCCESS;
    }
    else{
      if(Eleft>Eright){
	for(i=0;i<nvar;i++) min[i] = pleft[i];
      }
      else {
	for(i=0;i<nvar;i++) min[i] = pright[i];
      }
    }
  }

  printf("Needed more than %d cycles to find bracket on sphere\n",MAXCYCLE);
  printf("Leaving ivar=%d for now and continuing with next dimension..\n",ivar);
  return EXIT_FAILURE;
}

// -------------------------------------------------------------
void minimize_energy_onsphere(HILLSPATH *hills, double *min, double *center, double stepsize, double bracketsize, int ivar, FILE *fpout){

  double pright[NVARMAX], pleft[NVARMAX], point[NVARMAX], Epoint;
  double step, Eright, Eleft, Emiddle, fact;
  int  i, ihill, icycle, nvar;

  nvar=hills->nvar;
  Emiddle=0.0;
  for(i=0;i<nvar;i++) {
    pleft[i] = min[i] ;  pright[i] = min[i] ;  point[i] = min[i] ;
  }
  /* make right bracket point on sphere */
  pright[ivar] = pright[ivar] + bracketsize;
  //    print_point(fpout,hills,pright);
  fact=stepsize/radius(pright,center,nvar);
  for(i=0;i<nvar;i++) pright[i]= center[i]+(pright[i]-center[i])*fact;
  //  if(ivar==0) print_point(fpout,hills,pright);
  /* make left bracket point on sphere */
  pleft[ivar] = pleft[ivar] - bracketsize;
  //   print_point(fpout,hills,pleft);
  fact=stepsize/radius(pleft,center,nvar);
  for(i=0;i<nvar;i++) pleft[i]= center[i]+(pleft[i]-center[i])*fact;
  //  if(ivar==0) print_point(fpout,hills,pleft);
  Eleft=energy(pleft,hills);
  Eright=energy(pright,hills);
  Emiddle=energy(min,hills);

  icycle=0;
  while( ( fabs(pright[ivar]-pleft[ivar]) > TOL) && (icycle < 10000) ){
    //    printf("%d  %d    %lf %lf     %lf %lf     %lf %lf\n",ivar,icycle,pleft[ivar],Eleft,min[ivar],Emiddle,pright[ivar],Eright);
    if(Eleft>Eright){          /* pick new point between pleft and min */
      step=0.5*(pright[ivar] - min[ivar]);
      point[ivar] = pright[ivar] - step;
      fact=stepsize/radius(point,center,nvar);
      for(i=0;i<nvar;i++) point[i]= center[i]+(point[i]-center[i])*fact;
      //      if(ivar==0)print_point(fpout,hills,point);
      Epoint = energy(point,hills);
      if(Epoint<Emiddle){
	for(i=0;i<nvar;i++) pright[i]=point[i];
	Eright=Epoint;
      }
      else{
	for(i=0;i<nvar;i++){
	  pleft[i]=min[i];
	  min[i]=point[i];
	}
	Eleft=Emiddle;
        Emiddle=Epoint;
      } 
    }
    else{                      /* pick new point between pright and min */
      step=0.5*(min[ivar]-pleft[ivar]);
      point[ivar] = pleft[ivar] + step;
      fact=stepsize/radius(point,center,nvar);
      for(i=0;i<nvar;i++) point[i]= center[i]+(point[i]-center[i])*fact;
      //      if(ivar==0)print_point(fpout,hills,point);
      Epoint = energy(point,hills);
      if(Epoint<Emiddle){
	for(i=0;i<nvar;i++) pleft[i]=point[i];
	Eleft=Epoint;
      }
      else{
	for(i=0;i<nvar;i++){
	  pright[i]=min[i];
	  min[i]=point[i];
	}
	Eright=Emiddle;
	Emiddle=Epoint;
      }
    }
    //    printf("cycle %d: pright-pleft=%e\n",icycle,pright[ivar]-pleft[ivar]);
    icycle++;
  }
  if(icycle > 10000) printf("WARNING minimization stopped after %d cycles. (pright[ivar]-pleft[ivar]=%lg)\n",icycle,pright[ivar]-pleft[ivar]);
  //  fprintf(fpout,"\n");
  //  printf("> %d %lf",icycle,min[0]);
  //  for(i=1;i<nvar;i++) printf(" %lf",min[i]);
  //  printf("  %lg\n",energy(min,hills));
}

// -------------------------------------------------------------
int check_trace_direction(double *prevmin, double *min, double *point, int nvar, int npoints,int *tstepdirection){

  int i;

  if( radius(prevmin,point,nvar) < DISTTOL ) {
    printf("Distance between point %d and point %d is very small\n",npoints+1, npoints-1);
    printf("Seems like we are going back\n");
    if(*tstepdirection){
      printf("Already tried resetting point to opposite direction\n");
      printf("Program stopped\n");
      return EXIT_FAILURE;
    }
    else{
      printf("Resetting point in opposite direction\n");
      *tstepdirection=1;
      for(i=0;i<nvar;i++) point[i]=2*min[i]-prevmin[i];
    }
  }
  return EXIT_SUCCESS;
}

// -------------------------------------------------------------
void crossproduct(double *vec1, double *vec2, double *vec3){

  double dummy[3];
  int i;

  for(i=0;i<3;i++) dummy[i]=vec3[i];
  dummy[0]=vec1[1]*vec2[2]-vec1[2]*vec2[1];
  dummy[1]=vec1[2]*vec2[0]-vec1[0]*vec2[2];
  dummy[2]=vec1[0]*vec2[1]-vec1[1]*vec2[0];
  for(i=0;i<3;i++) vec3[i]=dummy[i];
}



// -------------------------------------------------------------
void get_perpendicular_higher(double *vec, double *perpvec, int ivar, int nvar, double length){

  double origin[nvar];
  double perplength;
  int i,jvar;

  /* make perpendicular vector by setting all elements to zero, except
     element ivar (set to 1), and solving for element jvar=ivar+1 so that
     the dot product between vec and perpvec equals zero */
  for(i=0;i<nvar;i++) perpvec[i] = 0.0;
  for(i=0;i<nvar;i++) origin[i] = 0.0;
  jvar=ivar+1;
  if(jvar == nvar) jvar=0;
  if(fabs(vec[jvar])>0.00001){
    perpvec[ivar] = 1.0;
    perpvec[jvar] = -1.0*vec[ivar]/vec[jvar];
  }else{
    if(fabs(vec[ivar])>0.00001){
      perpvec[jvar] = 1.0;
      perpvec[ivar] = -1.0*vec[jvar]/vec[ivar];
    }
    else{
      perpvec[jvar] = 1.0;
      perpvec[ivar] = 0.0;
    }
  }
  perplength = radius(perpvec,origin,nvar);
  for(i=0;i<nvar;i++) perpvec[i]=perpvec[i]*length/perplength;
  return;
}

// -------------------------------------------------------------
void get_perpendicular_3Dvector(double *vec, double *perpvec, double length){

  double origin[3]={0.0, 0.0, 0.0};
  double perplength;
  int i;

  /* make perpendicular vector by taking the crossproduct of vec and a eenheidsvektor e1 */
  perpvec[0] = 0.0;
  perpvec[1] = vec[2];
  perpvec[2] = -vec[1];
  /* perhaps vec was almost parallel to this eenheidsvector e1 */
  perplength = radius(perpvec,origin,3);
  if(perplength < 0.05*length){
    perpvec[0] = -vec[2];
    perpvec[1] = 0.0;
    perpvec[2] = vec[0];
    perplength = radius(perpvec,origin,3);
  }
  for(i=0;i<3;i++) perpvec[i]=perpvec[i]*length/perplength;
}


// -------------------------------------------------------------
int bracket_minimim_online(HILLSPATH *hills, double *min, double *pvec, double bracketsize, FILE *fpout, int icolor){

  double pointright[NVARMAX], pointleft[NVARMAX];
  double Eright, Eleft, Emiddle;
  int ivar,i;


  //  printf("============= Bracketing minimum ===============================\n");
  //  printf("    left bracket   Energy           minimum      Energy       right bracket  Energy \n");
  for(ivar=0;ivar<hills->nvar;ivar++){
    pointleft[ivar]  = min[ivar] + bracketsize*pvec[ivar];
    pointright[ivar] = min[ivar] - bracketsize*pvec[ivar];
  }
  //	print_point(fpout,hills,pointleft);
  //	print_point(fpout,hills,pointright);
 /* TCL/TK SCRIPT OUTPUT */
#if PRINT_TCL_BRACKET==1
  while(icolor>=10) icolor-=10;
  fprintf(fpout,"plot line  ");
  for(i=0;i<hills->nvar;i++) fprintf(fpout,"%lf ",min[i]);
  for(i=0;i<hills->nvar;i++) fprintf(fpout,"%lf ",pointleft[i]);
  fprintf(fpout,"%s append\n",color[icolor]);
  fprintf(fpout,"plot line  ");
  for(i=0;i<hills->nvar;i++) fprintf(fpout,"%lf ",min[i]);
  for(i=0;i<hills->nvar;i++) fprintf(fpout,"%lf ",pointright[i]);
  fprintf(fpout,"%s append\n",color[icolor]);
#endif
 /* END OF TCL/TK SCRIPT OUTPUT */
  Eleft=energy(pointleft,hills);
  Eright=energy(pointright,hills);
  Emiddle=energy(min,hills);
  // printf("  %lf %lf %lf   %lf %lf %lf   %lf %lf %lf\n",pointleft[0],pointleft[1],Eleft,min[0],min[1],Emiddle,pointright[0],pointright[1],Eright);
  if((Eleft>Emiddle)&&(Eright>Emiddle)){
    printf("Bracketing failed.\nAdjust bracketsize.\nProgram stopped.\n");
    return EXIT_FAILURE;
  }
  if(Eleft>Emiddle){
    for(ivar=0;ivar<hills->nvar;ivar++) min[ivar] = pointleft[ivar];
  }
  else{
    for(ivar=0;ivar<hills->nvar;ivar++) min[ivar] = pointright[ivar];
  }
  return EXIT_SUCCESS;
}


// -------------------------------------------------------------
void minimize_energy_online(HILLSPATH *hills, double *min, double *pvec, double bracketsize, FILE *fpout, int icolor){

  double pright[NVARMAX], pleft[NVARMAX], point[NVARMAX], Epoint;
  double step, Eright, Eleft, Emiddle, Eprevious, fact;
  int  i, icycle, nvar;

  nvar=hills->nvar;
  Emiddle=0.0;
  Eprevious=Emiddle-2.0*ETOL;
  Eprevious=Emiddle;
  for(i=0;i<nvar;i++){
    pleft[i]  = min[i] + bracketsize*pvec[i];
    pright[i] = min[i] - bracketsize*pvec[i];
  }
  Eleft=energy(pleft,hills);
  Eright=energy(pright,hills);
  Emiddle=energy(min,hills);

  icycle=0;
  //  printf("cycle %d: pright-pleft=%e\n",icycle,radius(pright,pleft,nvar));
  while( (radius(pleft,pright,nvar)  > TOL) && (icycle < 100) ){
    if(Eleft<Eright){          /* pick new point between pleft and min */
      for(i=0;i<nvar;i++) point[i]=0.5*(pright[i] + min[i]);
      Epoint = energy(point,hills);
      if(Epoint<Emiddle){
	for(i=0;i<nvar;i++) pright[i]=point[i];
	Eright=Epoint;
      }
      else{
	for(i=0;i<nvar;i++){
	  pleft[i]=min[i];
	  min[i]=point[i];
	}
	Eleft=Emiddle;
	Emiddle=Epoint;
      } 
    }
    else{                      /* pick new point between pright and min */
      for(i=0;i<nvar;i++) point[i]=0.5*(pleft[i] + min[i]);
      Epoint = energy(point,hills);
      if(Epoint<Emiddle){
	for(i=0;i<nvar;i++) pleft[i]=point[i];
        Eleft=Epoint;
      }
      else{
	for(i=0;i<nvar;i++){
	  pright[i]=min[i];
	  min[i]=point[i];
	}
	Eright=Emiddle;
	Emiddle=Epoint;
      }
    }
    icycle++;
    //    printf("cycle %d: pright-pleft=%e\n",icycle,radius(pright,pleft,nvar));
    //    print_point(fpout,hills,point);
  }
  if(icycle > 10000) printf("WARNING minimization stopped after %d cycles. (distance(pright-pleft)=%lg)\n",icycle,radius(pright,pleft,nvar));
  //  fprintf(fpout,"\n");
#if PRINT_TCL_MINIMIZE==1
  /* TCL/TK SCRIPT OUTPUT */
  fprintf(fpout,"plot sphere  ");
  for(i=0;i<hills->nvar;i++) fprintf(fpout,"%lf ",point[i]);
  while(icolor>=10) icolor-=10;
  fprintf(fpout,"0.01 %s append\n",color[icolor]);
  /* END OF TCL/TK SCRIPT OUTPUT */
#endif
  //  printf("> %d %lf",icycle,min[0]);
  //  for(i=1;i<nvar;i++) printf(" %lf",min[i]);
  //  printf("  %lg\n",energy(min,hills));
}

 // -------------------------------------------------------------
void get_coursegrid_vector(double **points, int nvar, int istep, int ivar, int ipoint, int nsteps,int npoints, double length, double *vec, double *vec2){

  /* vec is the vector of course grid points[i+1]-points[i]  */
  /* vec2 an average vector of points[i+1]-points[i] and points[i]-points[i-1] */

  double vec1[NVARMAX],vec3[NVARMAX];

  for(ivar=0;ivar<nvar;ivar++) {
    vec[ivar]=(points[ipoint+1][ivar]-points[ipoint][ivar])/length;
    if( (istep <= nsteps/2) && (ipoint!=0) ){
      vec1[ivar]=(points[ipoint][ivar]-points[ipoint-1][ivar]) /
	radius(points[ipoint],points[ipoint-1],nvar);
      vec2[ivar]=((nsteps/2-istep)*vec1[ivar]+(nsteps/2+istep)*vec[ivar])/nsteps;
    }
    else if( (istep > nsteps/2) && (ipoint!=(npoints-2)) ){
      vec3[ivar]=(points[ipoint+2][ivar]-points[ipoint+1][ivar]) /
	radius(points[ipoint+2],points[ipoint+1],nvar);
      vec2[ivar]=((3*nsteps/2-istep)*vec[ivar]+(istep-nsteps/2)*vec3[ivar])/nsteps;
    }
    else{
      vec2[ivar]=vec[ivar];
    }
  }
}

 // -------------------------------------------------------------
double integrate_perp2irc(HILLSPATH *hills, double *min, double *vec, FILE *fpout){
  int i;
  double free_energy, stepsize, point[2], fact;

  stepsize=0.005;
  for(i=-50;i<51;i++){
    point[0] = min[0] + i*stepsize*vec[1];
    point[1] = min[1] + i*stepsize*vec[0];
    if(i==0){
      fact=0.0;
    }
    else{
      fact=i/abs(i);
    }
    fprintf(fpout,"  %d  %lf  ",i,fact*radius(point,min,hills->nvar));
    print_point(fpout,hills,point);
  }
  return free_energy;
}

// -------------------------------------------------------------

int get_irc_value(int nvar,int totpoints, double **irc_traj, double *point, double *radius){

  /* THIS PROBABLY GOES WRONG FOR COLL.VARS OF DIFFERENT DIMENSION !!! */
  /* YOU MIGHT WANT TO SCALE THE VARIABLES FIRST TO SOME SORT OF GENERALIZED COORDS */

  double dxmin,dx,diff2;
  int i,ipoint,irc_val;

  dxmin=999999999999.;
  for(ipoint=0;ipoint<totpoints;ipoint++){
    dx=0.0;
    for(i=0;i<nvar;i++){
      diff2=irc_traj[ipoint][i]-point[i];
      diff2*=diff2;
      dx+=diff2;
    }
    if(dx<dxmin){
      dxmin=dx;
      irc_val=ipoint;
    }
  }
  *radius=sqrt(dxmin);
  return irc_val;
}

// -------------------------------------------------------------
  /* Perform Monte Carlo simulation to compute an effective umbrella */
void montecarlo(HILLSPATH *hills, double *startpoint, int maxmcmoves, double mcmaxetot, double drmax, double dxmax, int totpoints, double **irc_traj,double *E_traj, double *DE_traj, FILE *fptcl){

  int i,ntrail, naccpt, nvar,tmoremoves,irc_val,ipoint;
  double vnew,vold,deltv,deltvb, ratio;
  const double beta=1052.582355;    /* 4184*627.509556/300.0/8.314472  */
  double oldp[NVARMAX], newp[NVARMAX], *histo, dx, pot, p_irc;
  FILE *fp;

  dseed=428769461.0;
  printf("\n === Entering Monte Carlo simulation ===\n");
  if( (fp = fopen(MCOUTPUTFILE,"w")) == NULL){      
    printf("Cannot open %s \nStopped\n",MCOUTPUTFILE);
    return;
  }
  nvar=hills->nvar;
  ntrail=0;
  naccpt=0;
  tmoremoves=1;
  histo = (double *) malloc(totpoints*sizeof(double));
  for(ipoint=0;ipoint<totpoints;ipoint++) histo[ipoint]=0;  
  for(i=0;i<nvar;i++) oldp[i]=startpoint[i];
  irc_val=get_irc_value(nvar,totpoints,irc_traj,oldp,&dx);
  vold=E_traj[irc_val]-energy(oldp,hills);
  while(tmoremoves){
    for(i=0;i<nvar;i++){
      newp[i]= oldp[i] + (2.0*myrandom(0)-1.0)*drmax;
    }
    irc_val=get_irc_value(nvar,totpoints,irc_traj,newp,&dx);
    pot=energy(newp,hills);
    vnew=E_traj[irc_val]-pot;
    //    vnew=-energy(newp,hills);
    deltv = vnew - vold;
    deltvb = beta * deltv;
    //    if( (deltvb < 75.0) && (dx < dxmax) ){   /* only sample neigborhood of IRC */
    if( (deltvb < 75.0) && (pot > mcmaxetot) ){  /* only sample significant hill population */
      if( deltvb<=0.0 ){
	naccpt++;
	histo[irc_val]++;
#if PRINT_MC_MOVES_TCL==1
        fprintf(fptcl,"plot sphere  ");
        for(i=0;i<nvar;i++) fprintf(fptcl,"%lf ",newp[i]);
        fprintf(fptcl,"0.01 green append\n");
/* 	fprintf(fptcl,"plot line  "); */
/* 	for(i=0;i<nvar;i++) fprintf(fptcl,"%lf ",oldp[i]); */
/* 	for(i=0;i<nvar;i++) fprintf(fptcl,"%lf ",newp[i]); */
/*         fprintf(fptcl,"green append\n"); */
#endif
	for(i=0;i<nvar;i++){
          oldp[i]=newp[i];
	}
	vold=vnew;
#if PRINT_MC_MOVES==1
      fprintf(fp,"%d %d %lf",naccpt,ntrail,drmax);
      for(i=0;i<nvar;i++){
	fprintf(fp," %lf",oldp[i]);
      }
      fprintf(fp," %lf\n",vold);
#endif
      }
      else if(exp(-deltvb)>myrandom(0)){
	naccpt++;
	histo[irc_val]++;
#if PRINT_MC_MOVES_TCL==1
        fprintf(fptcl,"plot sphere  ");
        for(i=0;i<nvar;i++) fprintf(fptcl,"%lf ",newp[i]);
        fprintf(fptcl,"0.01 blue append\n");
/* 	fprintf(fptcl,"plot line  "); */
/* 	for(i=0;i<nvar;i++) fprintf(fptcl,"%lf ",oldp[i]); */
/* 	for(i=0;i<nvar;i++) fprintf(fptcl,"%lf ",newp[i]); */
/*         fprintf(fptcl,"blue append\n"); */
#endif
	for(i=0;i<nvar;i++){
          oldp[i]=newp[i];
	}
	vold=vnew;
#if PRINT_MC_MOVES==1
	fprintf(fp,"%d %d %lf",naccpt,ntrail,drmax);
	for(i=0;i<nvar;i++){
	  fprintf(fp," %lf",oldp[i]);
	}
	fprintf(fp," %lf\n",vold);
#endif
      }
    }
    ntrail++;
    /* tune stepsize */
#if TUNE_MC_STEPSIZE==1
    if(fmod(ntrail,NADJUST)==0){
      ratio = (double)naccpt / (double)ntrail;
      if(ratio > 0.5){
	drmax*=1.05;
      }
      else{
	drmax*=0.95;
      }
    }
#endif
    if(naccpt>maxmcmoves) tmoremoves=0;
  }

  fclose(fp);
  /* write histogram to file and compute dE = kBT ln(P(IRC) */
  if( (fp = fopen(MCHISTOOUT,"w")) == NULL){      
    printf("Cannot open %s \nStopped\n",MCHISTOOUT);
    return;
  }
  for(ipoint=0;ipoint<totpoints;ipoint++){
    fprintf(fp,"%d ",ipoint);
    for(i=0;i<nvar;i++) fprintf(fp," %lf",irc_traj[ipoint][i]);
    p_irc = (double)histo[ipoint]/(double)naccpt;
    DE_traj[ipoint]= -KBT*log(p_irc);
    fprintf(fp," %lf  %lf\n",p_irc,DE_traj[ipoint]);
  }
  free(histo);
  fclose(fp);
  printf("\n === Finished Monte Carlo simulation ===\n");
}




// -------------------------------------------------------------

int main(int argc, char *argv[]){

  FILE *fpout, *fpout2, *fpout3, *fpout4;
  HILLSPATH hills;
  double **points;   /* [NPOINTSMAX][NVARMAX]  stores first and second minumum and points in between */
  double point[NVARMAX];
  double bracketsize, stepsize, Eprevious, Epoint;
  int status, i, ivar, npoints, tedgeconst, noffset, icycle, tmc, maxmcmoves;
  int tstepdirection;   /* =0 for step in direction of other minimum and
                           =1 for step in same direction of previous step  */
  int istep, nsteps, ipoint, icount, jcount, totpoints, nside, numberofsteps;
  int tfirst, itotal, tminopt, ttrace;
  int inputtype, nprint_fesmovie;
  double vec[NVARMAX],avervec[NVARMAX],perpvec[NVARMAX];
  double free_energy, length, cplane, fact, fact2, dx, mcminetot;
  double  **irc_traj, gausswidth[2], norm, weight, *E_traj, *DE_traj;
  double force[NVARMAX], **temp, nedgepnt, *Etemp, maxstep, dxmax;
  char colvarfile[32],hillfile[32];

  /* print helpfile */
  if(argc>1){
    print_help();
    return EXIT_FAILURE;
  }


  /* set defaults */
  inputtype=0;              /* CPMD=0, PLUMED=1, POTENTIAL=2  */
  tminopt=1;                /* if 0 then the minima are not further optimized */
  ttrace=1;                 /* if 0 then only the minima are localized */
  tmc=0;                    /* if 1 then a monte carlo simulation computes an effective */
  /*                        umbrella potential, integrating out DOF perpendiculat to IRC*/
  hills.cutgauss = 100000000.0;
  hills.tspherical = 0;       /* if 1 then spherical hills are being used */
  hills.tnewcode = 0;         /* if 1 then new format of parvar_mtd input file is used */
  hills.tsym = 0;             /* if 1 then symmetry info is expected in colvar files */
  hills.isym = -1;            /* if not -1 then only hills with isym are read */
  hills.nvar = -1;
  hills.nx1=-1;               /* grid size of input potential */
  hills.nx2=-1;               /* grid size of input potential */
  maxstep = 0.1;              /* max stepsize of Metropolis Monte Carlo move */
  dxmax = 999999999.;         /* max distance of MC trail moves from IRC */
  mcminetot = -99999999.;     /* minimum total energy trail moves can make to prevent */
        /* running away of MC simulation when hills are almost zero (typically in TS)*/
  bracketsize = 0.2;
  stepsize = 0.2;
  tedgeconst=0;
  nsteps=10;                 /* number of points between interval */
  nside=10;                  /* number of points on each side, outside of the minima */
  nprint_fesmovie=-1;        /* if >0 print evolution of FES with this interval of HILLs */
  gausswidth[0] = -999.;     /* factor for gaussian smoothing of the irc points */
  gausswidth[1] = -999.;     /* factor for gaussian smoothing of the energy */
  strcpy(colvarfile,COLVARFILE);
  strcpy(hillfile,DCOLVARFILE);
  points = (double **) malloc(NPOINTSMAX*sizeof(double *));
  for(i=0;i<NPOINTSMAX;i++) points[i] = (double *) malloc(NVARMAX*sizeof(double));


  /* read input */
  status = read_input(&hills, points, &bracketsize, &stepsize, &nsteps, gausswidth, &nside, &tedgeconst, &tminopt, &ttrace, &tmc, &maxmcmoves, &maxstep, &dxmax, &mcminetot,colvarfile,hillfile,&inputtype,&nprint_fesmovie);
  if(status) return EXIT_FAILURE;

  /* read hills from colvar_val and dcolvar_val or PLUMED HILLS file or potential on a grid */
  if(inputtype==0){
    energy = &energy_meta;    /* set energy pointer to metadynamics potential function */
    status = read_hills_from_CPMDfile(&hills,colvarfile,hillfile);
  }else if(inputtype==1){
    energy = &energy_meta;   /* set energy pointer to metadynamics potential function */
    status = read_hills_from_PLUMEDfile(&hills,hillfile);
  }else if(inputtype==2){
    energy = &energy_pot;    /* set energy pointer to 2Dpotential interpolation function */
    status = read_2Dpotential(&hills,hillfile);
    get_2Dextpot_derivatives(&hills);
  }
  if(status) return EXIT_FAILURE;

  /* open output files */
  if(ttrace==1){
    if( (fpout = fopen(OUTPUTFILE,"w")) == NULL){       /* irc on course grid      */
      printf("Cannot open %s \nStopped\n",OUTPUTFILE);
      return EXIT_FAILURE;
    }
    if( (fpout2 = fopen(OUTPUTFILE2,"w")) == NULL){     /* bracket and optim. info */
      printf("Cannot open %s \nStopped\n",OUTPUTFILE2);
      return EXIT_FAILURE;
    }
    if( (fpout3 = fopen(OUTPUTFILE3,"w")) == NULL){     /* irc on small grid       */
      printf("Cannot open %s \nStopped\n",OUTPUTFILE3);
      return EXIT_FAILURE;
    }
    if( (fpout4 = fopen(OUTPUTFILE4,"w")) == NULL){     /* integral perp2irc info  */
      printf("Cannot open %s \nStopped\n",OUTPUTFILE4);
      return EXIT_FAILURE;
    }
  }

  
  /* minima */
  printf("\n  Energy in guess for minimum1 [%lf",points[0][0]);
  for(i=1;i<hills.nvar;i++) printf(",%lf",points[0][i]);
  printf("]= %lg\n",energy(points[0],&hills));
  if(tminopt){     /* optimize minima */
    status = bracket_minimum(&hills, points[0], bracketsize);
    if(status) return EXIT_FAILURE;

    printf("bracket minimum1 ok\n");
    minimize_energy(&hills, points[0], bracketsize);
    printf("  Minimum1 found. Energy for minimum1 [%lf",points[0][0]);
    for(i=1;i<hills.nvar;i++) printf(",%lf",points[0][i]);
    printf("]= %lg\n",energy(points[0],&hills));
  }
  
  printf("\n  Energy in guess for minimum2 [%lf",points[1][0]);
  for(i=1;i<hills.nvar;i++) printf(",%lf",points[1][i]);
  printf("]= %lg\n",energy(points[1],&hills));
  if(tminopt){     /* optimize minima */
    status = bracket_minimum(&hills, points[1], bracketsize);
    if(status) return EXIT_FAILURE;

    printf("bracket minimum2 ok\n");
    minimize_energy(&hills, points[1], bracketsize);
    printf("  Minimum2 found. Energy for minimum2 [%lf",points[1][0]);
    for(i=1;i<hills.nvar;i++) printf(",%lf",points[1][i]);
    printf("]= %lg\n\n",energy(points[1],&hills));
  }
  
  //  print_point(fpout,&hills,points[1]);
  //  print_point(fpout,&hills,points[0]);

  if(ttrace==0){       /* no trace, only minima localization */
    printf("Program finished succesfully\n");
    return(EXIT_SUCCESS);
  }

  /* COURSE GRID COURSE GRID COURSE GRID COURSE GRID COURSE GRID COURSE GRID */
  /* trace IRC on a course grid first (only if the number of variable is more than one) */  

  npoints=0;
  if(hills.nvar>1){
    printf("  Tracing course path\n");
    for(i=0;i<hills.nvar;i++) point[i]=points[0][i];  /* startpoint=first element */

    while( radius(points[npoints+1],point,hills.nvar) > stepsize ){
      //  while(npoints<2){
      /* pick new point on the line from previous point and min2 */
      tstepdirection=0;
      for(i=0;i<hills.nvar;i++) {
	point[i]+=( points[npoints+1][i] - points[npoints][i]) * 
	  stepsize/radius(points[npoints+1],points[npoints],hills.nvar);
      }

      //    fprintf(fpout,"\n");
      //      print_point(fpout,&hills,point);

      Eprevious=0.0;
      Epoint=energy(point,&hills);
      while( (fabs(Epoint-Eprevious)) > ETOL ){
	for(ivar=0;ivar<hills.nvar;ivar++){
	  status = bracket_minimum_onsphere(&hills, point, points[npoints], 
					    stepsize, bracketsize, ivar, fpout2);
	  //	if(status) return EXIT_FAILURE;
	  if(status) continue;   // if status: did not find bracket, skipping to next var
	  minimize_energy_onsphere(&hills, point, points[npoints], stepsize, bracketsize, ivar, fpout2);
	  if(npoints > 0) status = check_trace_direction(points[npoints-1], 
			       points[npoints], point, hills.nvar,npoints,&tstepdirection);
	  if(status) return EXIT_FAILURE;
	}
	Eprevious=Epoint;
	Epoint=energy(point,&hills);
      }
      //      print_point(fpout,&hills,point);
      npoints++;
      printf("npoint= %d %lf",npoints,point[0]);
      for(i=1;i<hills.nvar;i++) printf(" %lf",point[i]);
      printf("  %lg\n",energy(point,&hills));
      /* found point becomes new startpoint and minimum2 is shifted one up*/    
      for(i=0;i<hills.nvar;i++){
	if(npoints >= NPOINTSMAX){
	  printf("Number of points is more than %d\n",NPOINTSMAX);
	  printf("Increase NPOINTSMAX\nProgram stopped;");
	  return EXIT_FAILURE;
	}
	points[npoints+1][i]=points[npoints][i];
	points[npoints][i]=point[i];
      }
    }
    npoints+=2;
  /* print course grid point to file and stdout */
    printf("EUREKA!\nConnected the two minima with %d points\n",npoints-2);
    printf("(these points are also written to %s)\n",OUTPUTFILE);
    for(i=0;i<npoints;i++){
      printf(" %d  ",i);
      print_point(stdout,&hills,points[i]);
      print_point(fpout,&hills,points[i]);
    }
  }

  /* print vector data to trace_out.tcl for gOpenmol */
  fprintf(fpout4,"3 200\n");
  fprintf(fpout4,"%d  %d  %d\n",npoints-1,npoints-1,npoints-1);
  fprintf(fpout4,"0.00000 1.50000\n");
  fprintf(fpout4,"0.00000 1.50000\n");
  fprintf(fpout4,"1.75000 3.25000\n");
  for(i=0;i<npoints-1;i++){
    fprintf(fpout4,"%lf %lf %lf %lf %lf %lf\n",points[i][0],points[i][1],points[i][2],points[i+1][0]-points[i][0],points[i+1][1]-points[i][1],points[i+1][2]-points[i][2]);
  }
  fflush(fpout4);


  /* FINE GRID FINE GRID FINE GRID FINE GRID FINE GRID FINE GRID */

  /* trace IRC between the previously found points by         */
  /* trying to locate the minimum Energy in the hyperplane    */
  /* perpendicular to the vecor from points[i] to point[i+1]  */
  /* Plane: v1 x1 + v2 x2 + .. + Const =0, {v1..vn}=vector    */


  printf("\nEntering fine grid IRC search\nBracketsize=%lf\n",bracketsize);
  fprintf(fpout,"\n");
  icount=0;
  itotal=0;
  tfirst=1;
  irc_traj = (double **) malloc((2*nside+npoints*nsteps)*sizeof(double *));
  for(i=0;i<(2*nside+npoints*nsteps);i++) irc_traj[i]=(double *) malloc(hills.nvar*sizeof(double));
  printf("Allocating memory for %d points.\n",2*nside+npoints*nsteps);
  dx=radius(points[1],points[0],hills.nvar)/nsteps;
  for(ipoint=0;ipoint<npoints-1;ipoint++){
    length=radius(points[ipoint+1],points[ipoint],hills.nvar);
    istep=0;
    if(ipoint == 0) istep=-nside;    /* start with nside points before minimum */
    numberofsteps=nsteps;            /* is increased to get nside points extra at end */
    while( istep<numberofsteps){
      /* vec is the vector of course grid points[i+1]-points[i]  */
      /* avervec an average vector of points[i+1]-points[i] and points[i]-points[i-1] */
      get_coursegrid_vector(points,hills.nvar,istep,ivar,ipoint,nsteps,npoints,length,vec,avervec); 
      fact=dx*istep;
      switch(hills.nvar){
      case 1:                            /* for one dimension */
	irc_traj[icount][0]=points[0][0]+fact;
        break;
      case 2:                            /* for two dimensions */
        for(i=0;i<hills.nvar;i++) point[i]=points[ipoint][i] + fact*vec[i];
        perpvec[0] = avervec[1];
	perpvec[1] = -avervec[0];
 	status = bracket_minimim_online(&hills, point, 
					perpvec, bracketsize, fpout2,0);
	if(status) return EXIT_FAILURE;
	minimize_energy_online(&hills, point, perpvec, 
					bracketsize, fpout2,0);
        for(i=0;i<hills.nvar;i++) irc_traj[icount][i]=point[i];
	//	print_point(fpout3,&hills,point);
	//      free_energy=integrate_perp2irc(&hills, point, avervec, fpout4);
	break;
      case 34:                           /* for three dimensions */
        for(i=0;i<hills.nvar;i++) point[i]=points[ipoint][i] + fact*vec[i];
	get_perpendicular_3Dvector(avervec,perpvec,dx*(nsteps+1));
        Eprevious=0.0;
        Epoint=energy(point,&hills);
        icycle=0;
        while( (fabs(Epoint-Eprevious)> ETOL) && ( icycle < MAXCYCLE) ){
	  status = bracket_minimim_online(&hills, point,
					  perpvec, bracketsize, fpout2,icycle);
	  if(status) return EXIT_FAILURE;
	  minimize_energy_online(&hills, point, perpvec, 
				 bracketsize, fpout2, icycle);
	  crossproduct(avervec,perpvec,perpvec);
  	  Eprevious=Epoint;
	  Epoint=energy(point,&hills);
	  icycle++;
	}
        printf("optimized point in %d cycles:",icycle);
        print_point(stdout,&hills,point);
        for(i=0;i<hills.nvar;i++) irc_traj[icount][i]=point[i];
	break;
      default:                    /* for higher dimensions (and perhaps also lower) */
        for(i=0;i<hills.nvar;i++) point[i]=points[ipoint][i] + fact*vec[i];
        Eprevious=0.0;
        Epoint=energy(point,&hills);
        icycle=0;
        while( (fabs(Epoint-Eprevious)> ETOL) && ( icycle < MAXCYCLE) ){
	  for(i=0;i<hills.nvar;i++){
	    get_perpendicular_higher(avervec,perpvec,i,hills.nvar,dx*(nsteps+1));
	    status = bracket_minimim_online(&hills, point,
					    perpvec, bracketsize, fpout2,icycle);
	    if(status) return EXIT_FAILURE;
	    minimize_energy_online(&hills, point, perpvec, 
				   bracketsize, fpout2, icycle);
	    Eprevious=Epoint;
	    Epoint=energy(point,&hills);
	    icycle++;
	  }
	}
        printf("optimized point in %d cycles:",icycle);
        print_point(stdout,&hills,point);
        for(i=0;i<hills.nvar;i++) irc_traj[icount][i]=point[i];
/* 	printf("Higher dimensions than 2 not yet implemented\nProgram stopped\n"); */
/* 	return EXIT_FAILURE; */
      }
      if(tfirst && ((radius(irc_traj[icount],points[npoints-1],hills.nvar)<dx) || 
		(icount == nside+(npoints-1)*nsteps-1) ) ){
	numberofsteps=istep+nside+2;
        tfirst=0;
      }
      icount++;
      istep++;
    }
    printf("Optimized %d points between found course grid points %d and %d\n",icount-itotal,ipoint,ipoint+1);
    itotal=icount;
  }
  totpoints=icount;
  printf("Total number of points=%d\n",totpoints);


 /* TCL/TK SCRIPT OUTPUT */
#if PRINT_TCL_CLEANPOINTS==1
  for(ipoint=1;ipoint<totpoints;ipoint++){
    fprintf(fpout2,"plot line  ");
    for(i=0;i<hills.nvar;i++) fprintf(fpout2,"%lf ",irc_traj[ipoint-1][i]);
    for(i=0;i<hills.nvar;i++) fprintf(fpout2,"%lf ",irc_traj[ipoint][i]);
    fprintf(fpout2,"blue append\n");
  }
#endif
 /* END OF TCL/TK SCRIPT OUTPUT */


  /* smooth points  */

  //  nedgepnt=nside/2;               /* not smooting outher points "nedgepnt"  */
  nedgepnt=0;
  if(gausswidth[0]!=-999.){
    temp = (double **) malloc(totpoints*sizeof(double *));
    for(i=0;i<totpoints;i++) temp[i] =(double *)malloc(hills.nvar*sizeof(double));
    printf("Applying some data smoothing\n");
    gausswidth[0]*=gausswidth[0];
    for(icount=0;icount<nedgepnt;icount++){
      temp[icount]=irc_traj[icount];
    }
    for(icount=nedgepnt;icount<totpoints-nedgepnt;icount++){
      for(ivar=0;ivar<hills.nvar;ivar++){
	norm=0.0;
	temp[icount][ivar]=0.0;
	for(jcount=0;jcount<totpoints;jcount++){
	  weight = exp(-1.0*pow((icount-jcount),2)/gausswidth[0]);
	  norm+=weight;
	  temp[icount][ivar]+=irc_traj[jcount][ivar]*weight;
	}
	temp[icount][ivar]/=norm;
      }
    }
    for(icount=totpoints-nedgepnt;icount<totpoints;icount++){
      temp[icount]=irc_traj[icount];
    }
    free(irc_traj);
    irc_traj=temp;
    /* TCL/TK SCRIPT OUTPUT */
#if PRINT_TCL_SMOOTHPOINTS==1
    for(ipoint=1;ipoint<totpoints;ipoint++){
      fprintf(fpout2,"plot line  ");
      for(i=0;i<hills.nvar;i++) fprintf(fpout2,"%lf ",irc_traj[ipoint-1][i]);
      for(i=0;i<hills.nvar;i++) fprintf(fpout2,"%lf ",irc_traj[ipoint][i]);
      fprintf(fpout2,"red append\n");
    }
#endif
    /* END OF TCL/TK SCRIPT OUTPUT */
  }
  else{
    printf("No data smoothing applied\n");
  }
  printf("Total number of points=%d\n",totpoints);

  /* compute energy along IRC  */
  printf("Computing forces and writing stuff to file\n");
  E_traj = (double *) malloc(totpoints*sizeof(double ));
  if(tedgeconst){
    noffset=nside;
    printf("%d extra points added on both sides with constant energy (tedgeconst=true)\n",nside);
  }
  else{
    noffset=0;
    printf("%d extra points added on both sides with actual energy (tedgeconst=false)\n",nside);
  }
  for(ipoint=0;ipoint<noffset;ipoint++) 
    E_traj[ipoint]=energy(irc_traj[noffset],&hills);
  for(ipoint=noffset;ipoint<totpoints-noffset;ipoint++) 
    E_traj[ipoint]=energy(irc_traj[ipoint],&hills);
  for(ipoint=totpoints-noffset;ipoint<totpoints;ipoint++) 
    E_traj[ipoint]=energy(irc_traj[totpoints-noffset-1],&hills);


  /* Perform Monte Carlo simulation to compute an effective umbrella */
  DE_traj = (double *) malloc(totpoints*sizeof(double ));  
  for(ipoint=0;ipoint<totpoints;ipoint++) DE_traj[ipoint]=0.0;
  if(tmc == 1){
    if(tedgeconst==0){
      printf("\nWARNING YOU DID NOT SET TEDGECONST=1\nTHE MC SIMULATION CAN RUN AWAY\n");
    }
    ipoint=totpoints/2;
    montecarlo(&hills,irc_traj[ipoint],maxmcmoves,mcminetot,maxstep,dxmax,totpoints,irc_traj,E_traj,DE_traj,fpout2);
  }

  /* smooth energy  */
  if(gausswidth[1]!=-999.){
    gausswidth[1]*=gausswidth[1];
    printf("Applying some energy smoothing\n");
    Etemp = (double *) malloc(totpoints*sizeof(double ));
    nedgepnt=0;
    for(icount=nedgepnt;icount<totpoints-nedgepnt;icount++){
      norm=0.0;
      Etemp[icount]=0.0;
      for(jcount=0;jcount<totpoints;jcount++){
	weight = exp(-1.0*pow((icount-jcount),2)/gausswidth[1]);
	norm+=weight;
	Etemp[icount]+=(E_traj[jcount]-DE_traj[ipoint])*weight;
      }
      Etemp[icount]/=norm;
    }
    free(E_traj);
    E_traj=Etemp;
  }
  else{
    printf("No energy smoothing applied\n");
  }


  /* compute forces along IRC and print to OUTPUTFILE3  */
  ipoint=0;
  for(ivar=0;ivar<hills.nvar;ivar++) force[ivar]=0.0;  /* force in minimum1=0 */
  fprintf(fpout3,"  %d",ipoint+1);
  for(ivar=0;ivar<hills.nvar;ivar++) fprintf(fpout3,"   %lf",irc_traj[ipoint][ivar]);
  for(ivar=0;ivar<hills.nvar;ivar++) fprintf(fpout3,"   %lf",force[ivar]);
  fprintf(fpout3,"  %lf\n",E_traj[ipoint]);
  ipoint=0;
  fact2= 1.0 / radius(irc_traj[ipoint+1],irc_traj[ipoint],hills.nvar);
  fact2=fact2*fact2;
  fact2*=(E_traj[ipoint+1]-E_traj[ipoint]);
  for(ipoint=1;ipoint<totpoints-1;ipoint++){        /* force in points between minima */
    fact= 1.0 / radius(irc_traj[ipoint],irc_traj[ipoint+1],hills.nvar);
    fact=fact*fact;
    fact*=(E_traj[ipoint+1]-E_traj[ipoint]);
    for(ivar=0;ivar<hills.nvar;ivar++){
      force[ivar] =  fact2 * (irc_traj[ipoint][ivar]-irc_traj[ipoint-1][ivar]);
      force[ivar] += fact * (irc_traj[ipoint+1][ivar]-irc_traj[ipoint][ivar]);
      force[ivar] *= 0.5;
    }
    fact2 = fact;
    fprintf(fpout3,"  %d",ipoint+1);
    for(ivar=0;ivar<hills.nvar;ivar++) fprintf(fpout3,"   %lf",irc_traj[ipoint][ivar]);
    for(ivar=0;ivar<hills.nvar;ivar++) fprintf(fpout3,"   %lf",force[ivar]);
    fprintf(fpout3,"  %lf\n",E_traj[ipoint]);
  }
  ipoint=totpoints-1;
  for(ivar=0;ivar<hills.nvar;ivar++) force[ivar]=0.0;  /* force in minimum2=0 */
  fprintf(fpout3,"  %d",ipoint+1);
  for(ivar=0;ivar<hills.nvar;ivar++) fprintf(fpout3,"   %lf",irc_traj[ipoint][ivar]);
  for(ivar=0;ivar<hills.nvar;ivar++) fprintf(fpout3,"   %lf",force[ivar]);
  //  fprintf(fpout3,"  %lf\n",-1.0*E_traj[ipoint]);
  fprintf(fpout3,"  %lf\n",E_traj[ipoint]);


  if(nprint_fesmovie>0) print_movie(&hills,irc_traj,totpoints,nprint_fesmovie);

  printf("Program finished succesfully\n");
  fclose(fpout);
  fclose(fpout2);
  fclose(fpout3);
  fclose(fpout4);
  return(EXIT_SUCCESS);
}
