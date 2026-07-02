void rmatmult_3_at(double *dbl,int imin,
int imax,
int jmin,
int jmax,
int kmin,
int kmax,
int jp,
int kp,
 double *dbc,
 double *dbr,
 double *dcl,
 double *dcc,
 double *dcr,
 double *dfl,
 double *dfc,
 double *dfr,
 double *cbl,
 double *cbc,
 double *cbr,
 double *ccl,
 double *ccc,
 double *ccr,
 double *cfl,
 double *cfc,
 double *cfr,
 double *ubl,
 double *ubc,
 double *ubr,
 double *ucl,
 double *ucc,
 double *ucr,
 double *ufl,
 double *ufc,
 double *ufr, double *x, double *b )
{
   int i, ii, jj, kk ;

   /*
   int imin = *((int *)(rmatmult3_args[ 0 ]));
   int imax = *((int *)(rmatmult3_args[ 1 ]));
   int jmin = *((int *)(rmatmult3_args[ 2 ]));
   int jmax = *((int *)(rmatmult3_args[ 3 ]));
   int kmin = *((int *)(rmatmult3_args[ 4 ]));
   int kmax = *((int *)(rmatmult3_args[ 5 ]));
   int jp = *((int *)(rmatmult3_args[  6]));
   int kp = *((int *)(rmatmult3_args[  7]));
   */
   //printf("%d %d %d %d %d %d \n", imin, imax, jmin, jmax, kmin, kmax);
   /*
   const double *dbl = *((const double **)(rmatmult3_args[ 8 ]));
   const double *dbc = *((const double **)(rmatmult3_args[ 9 ]));
   const double *dbr = *((const double **)(rmatmult3_args[ 10 ]));
   const double *dcl = *((const double **)(rmatmult3_args[ 11 ]));
   const double *dcc = *((const double **)(rmatmult3_args[ 12 ]));
   const double *dcr = *((const double **)(rmatmult3_args[ 13 ]));
   const double *dfl = *((const  double **)(rmatmult3_args[ 14 ]));
   const double *dfc = *((const  double **)(rmatmult3_args[ 15 ]));
   const double *dfr = *((const  double **)(rmatmult3_args[ 16 ]));
   const double *cbl = *((const  double **)(rmatmult3_args[ 17 ]));
   const double *cbc = *((const  double **)(rmatmult3_args[ 18 ]));
   const double *cbr = *((const  double **)(rmatmult3_args[ 19 ]));
   const double *ccl = *((const  double **)(rmatmult3_args[ 20 ]));
   const double *ccc = *((const  double **)(rmatmult3_args[ 21 ]));
   const double *ccr = *((const  double **)(rmatmult3_args[ 22 ]));
   const double *cfl = *((const  double **)(rmatmult3_args[ 23 ]));
   const double *cfc = *((const  double **)(rmatmult3_args[ 24 ]));
   const double *cfr = *((const  double **)(rmatmult3_args[ 25 ]));
   const double *ubl = *((const  double **)(rmatmult3_args[ 26 ]));
   const double *ubc = *((const  double **)(rmatmult3_args[ 27 ]));
   const double *ubr = *((const  double **)(rmatmult3_args[ 28 ]));
   const double *ucl = *((const  double **)(rmatmult3_args[ 29 ]));
   const double *ucc = *((const  double **)(rmatmult3_args[ 30 ]));
   const double *ucr = *((const  double **)(rmatmult3_args[ 31 ]));
   const double *ufl = *((const  double **)(rmatmult3_args[ 32 ]));
   const double *ufc = *((const  double **)(rmatmult3_args[ 33 ]));
   const double *ufr = *((const  double **)(rmatmult3_args[ 34 ]));
   */
   
   double *xdbl = x - kp - jp - 1 ;
   double *xdbc = x - kp - jp     ;
   double *xdbr = x - kp - jp + 1 ;
   double *xdcl = x - kp      - 1 ;
   double *xdcc = x - kp          ;
   double *xdcr = x - kp      + 1 ;
   double *xdfl = x - kp + jp - 1 ;
   double *xdfc = x - kp + jp     ;
   double *xdfr = x - kp + jp + 1 ;
   double *xcbl = x      - jp - 1 ;
   double *xcbc = x      - jp     ;
   double *xcbr = x      - jp + 1 ;
   double *xccl = x           - 1 ;
   double *xccc = x               ;
   double *xccr = x           + 1 ;
   double *xcfl = x      + jp - 1 ;
   double *xcfc = x      + jp     ;
   double *xcfr = x      + jp + 1 ;
   double *xubl = x + kp - jp - 1 ;
   double *xubc = x + kp - jp     ;
   double *xubr = x + kp - jp + 1 ;
   double *xucl = x + kp      - 1 ;
   double *xucc = x + kp          ;
   double *xucr = x + kp      + 1 ;
   double *xufl = x + kp + jp - 1 ;
   double *xufc = x + kp + jp     ;
   double *xufr = x + kp + jp + 1 ;
   //double myflops = 0.0 ;
   
   
   //printf("calling the matmul \n");
   //i=6*jp + 7*kp;
   
   for ( kk = 2 ; kk < kmax ; kk++ ) {
       for ( jj = 2 ; jj < jmax ; jj++ ) {
           for ( ii = 2 ; ii < imax ; ii++ ) {
               i = ii + jj * jp + kk * kp ;
        
               b[i] = dbl[i] * xdbl[i] + dbc[i] * xdbc[i] + dbr[i] * xdbr[i] +
                   dcl[i] * xdcl[i] + dcc[i] * xdcc[i] + dcr[i] * xdcr[i] +
                   dfl[i] * xdfl[i] + dfc[i] * xdfc[i] + dfr[i] * xdfr[i] +
                   cbl[i] * xcbl[i] + cbc[i] * xcbc[i] + cbr[i] * xcbr[i] +
                   ccl[i] * xccl[i] + ccc[i] * xccc[i] + ccr[i] * xccr[i] +
                   cfl[i] * xcfl[i] + cfc[i] * xcfc[i] + cfr[i] * xcfr[i] +
                   ubl[i] * xubl[i] + ubc[i] * xubc[i] + ubr[i] * xubr[i] +
                   ucl[i] * xucl[i] + ucc[i] * xucc[i] + ucr[i] * xucr[i] +
                   ufl[i] * xufl[i] + ufc[i] * xufc[i] + ufr[i] * xufr[i] ;
	       
	       }
	       }
	       }
   
}



