#include "dpc_macros.h"
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <time.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

#define DPC_WARNING(STATUS,F1,F2) \
  printf("WARNING: %s ( in %s ) exited with status %d.\n",\
         #F1,#F2,STATUS);    
#define DPC_SING_WARNING printf("WARNING: matrix is singular in get_finalterm().\n");


int max_int(int a, int b){
  if(a > b){
    return(a);
  }else if(b > a){
    return(b);
  }else{
    return(a);
  }
}

double min_double(double a, double b){
  if(a < b){
    return(a);
  }else if(b < a){
    return(b);
  }else{
    return(a);
  }
}

void waitFor (unsigned int secs) {
    unsigned int retTime = time(0) + secs;   // Get finishing time.
    while (time(0) < retTime);               // Loop until it arrives.
}

/* Define Structures */ 
struct data {
  gsl_vector *k;
  gsl_vector *p;
  gsl_vector *prec;
  int *aInd;
  int n;
};

struct index {
  int start;
  int n;
};  


void brownian_precision(gsl_matrix * W){
  int i,j;
  int n_u = W->size1;
  for(i=0; i<n_u; i++){
    for(j=0; j<n_u; j++){
      if(i == j && (i == 0 || i == n_u -1)){
        gsl_matrix_set(W,i,j,1);
      }
      else if(i == j){
        gsl_matrix_set(W,i,j,2);
      }
      else if(i == j-1 || i == j + 1){
        gsl_matrix_set(W,i,j,-1);
      }
      else{
        gsl_matrix_set(W,i,j,0);
      }
      /* Print the matrix
      printf("%f ",gsl_matrix_get(W,i,j));
      if(j == *n_u - 1){ printf("\n"); }
      */
    }
  }
}


void data_creator(int *As, int *I, int *J, double *k, int N, 
                  int * NIs, int ni_max, int ntalldat,
                  double *p, double *prec, struct data DAT[N][ni_max],
                  int *dpc_status){

  struct index IND[N][ni_max];
  struct index tempind;
  int i,j,kk;
  int len = 0;
  int startind = 0;
  for(i=1; i < ntalldat; i++){
    len = len + 1;
    if( (I[i] != I[i-1] || J[i] != J[i-1]) || i == ntalldat - 1){ 
      tempind.n = len;
      tempind.start = startind;
      if(i < ntalldat - 1){ 
        IND[I[i-1]][J[i-1]] = tempind;
      }else{
        tempind.n = tempind.n + 1;
        IND[I[i]][J[i]] = tempind;
      }   
      startind = i;
      len = 0;
    }   
  }
  int cntr = 0;
  for(i=0; i < N; i++){
    for(j=0; j < NIs[i]; j++){
      gsl_vector * tempk = gsl_vector_alloc(IND[i][j].n);
      gsl_vector * tempp = gsl_vector_alloc(IND[i][j].n);
      gsl_vector * tempprec = gsl_vector_alloc(IND[i][j].n);
      if(!(tempk && tempp && tempprec)){
        ERR_MEM;
      }
      DAT[i][j].aInd = (int *)malloc(sizeof(int) * IND[i][j].n);
      if(!DAT[i][j].aInd){
        ERR_MEM;
      }
      /* Fill the new gsl array */
      cntr = 0;
      for(kk=IND[i][j].start; kk < IND[i][j].start + IND[i][j].n; kk++){
        gsl_vector_set(tempk,cntr,k[kk]);
        gsl_vector_set(tempp,cntr,p[kk]);
        gsl_vector_set(tempprec,cntr,prec[kk]);
        DAT[i][j].aInd[cntr] = As[kk];
        cntr = cntr + 1;
      }

      DAT[i][j].k = tempk;
      DAT[i][j].p = tempp;
      DAT[i][j].prec = tempprec;
      DAT[i][j].n = IND[i][j].n;
    }
  }
}

double get_C(int N, int *NIs, int ni_max, struct data DAT[N][ni_max]){
  int i,j,kk;
  double theC = 0;
  for(i=0; i < N; i++){
    for(j=0; j < NIs[i]; j++){
      for(kk=0; kk < DAT[i][j].n; kk++){
        // Obtain P^T Omega P for each dataset
        theC = theC - gsl_vector_get(DAT[i][j].prec,kk) *
                      gsl_vector_get(DAT[i][j].p,kk) *
                      gsl_vector_get(DAT[i][j].p,kk);
        // Obtain the log determinant of Omega
        theC = theC + log(gsl_vector_get(DAT[i][j].prec,kk));
      }
    }
  }
  theC = .5 * theC;
  return(theC);
}



double sum2(gsl_vector * x){
  int i;
  double out = 0.0;
  for(i=0; i < x->size; i++){
    out = out + gsl_vector_get(x,i) * gsl_vector_get(x,i);
  }
  return(out);
}


double rtruncnorm(const gsl_rng *r, double a, double b, double mu, double sigma){
  double val = 0.0;
  while(1){
    val = mu + gsl_ran_gaussian_ziggurat(r,sigma);
    if(val > a && val < b){
      return(val);
    }
  }
}

double dtruncnorm(double x, double a, double b, double mu, double sigma){
  double val = 0;
  val = gsl_ran_gaussian_pdf(x-mu,sigma);
  val = val / (gsl_cdf_gaussian_P(b-mu,sigma) - gsl_cdf_gaussian_P(a-mu,sigma));
  return(val);
}

/* NOTE: this function will modify A !!!*/ 
double get_finalterm(gsl_matrix *A, gsl_vector *x,int *status,int print_warn){
  #define ERR_MEM_GFT printf("ERROR: Not enough memory. File: %s (line %d)\n.",__FILE__,__LINE__); \
  *status = 2; \
  return(-1.0);

  double det,dval,out;
  int signum;
  gsl_permutation *p = gsl_permutation_alloc(A->size1);
  if(!p){ ERR_MEM_GFT }
  *status = gsl_linalg_LU_decomp(A,p,&signum);
  if(*status){
    if(print_warn){
      DPC_WARNING(*status,gsl_linalg_LU_decomp,get_finalterm());
    }
    return(-1.0);
  }
/*
  if(*status){
    printf("WARNING: gsl_linalg_LU_decomp (in get_finalterm()) \
           exited with status %d.\n",*status);
    return(-1.0);
  }
*/
  det = gsl_linalg_LU_lndet(A);
// printf("det=%f, ",det);
  gsl_vector * b = gsl_vector_alloc(x->size);
  if(!b){ ERR_MEM_GFT }
  *status = gsl_linalg_LU_solve(A,p,x,b); /* obtain C^{-1} x and store in b */
  if(*status == GSL_EDOM){
//    printf("WARNING: matrix is singular in get_finalterm().\n");
    if(print_warn){
      DPC_SING_WARNING
    }
    return(-1.0);
  }else if(*status){
    if(print_warn){
      DPC_WARNING(*status,gsl_linalg_LU_solve,get_finalterm())
    }
    return(-1.0);
  }
  *status = gsl_blas_ddot(x,b,&dval);
  if(*status){
     if(print_warn){
       DPC_WARNING(*status,gsl_blas_ddot,get_finalterm());
     }
     return(-1.0);
  }
  gsl_permutation_free(p);
  gsl_vector_free(b);
  out = .5*(dval - det);
// printf("dval=%f\n",dval);
  return out;
}

/* assumes the locs and aInd are based on sorted values  ORIGINAL FUNCTION 
void getAK(gsl_matrix *AK, int n_klocs, int n_u, gsl_matrix *K_sigma, int *aInd, int n_a){
  int i,j;
  int cntr = 0;
  int cntr_a = 0;
  for(i=0; i < n_klocs; i++){
    if(i == aInd[cntr_a]){
      for(j=0; j < n_u; j++){ 
        gsl_matrix_set(AK,cntr,j,gsl_matrix_get(K_sigma,i,j));
      } 
      cntr = cntr + 1;
      cntr_a = cntr_a + 1;
      if(cntr_a == n_a){
        return;
      }
    }
  }
}
*/ 

void getAK(gsl_matrix *AK, gsl_matrix *K_sigma, int *aInd){
  int i,j;
  int cntr = 0;
  int n_a = AK->size1;
  int n_u = K_sigma->size2;
//  int n_klocs = K_sigma->size1;
  /*
  printf("n_klocs=%i, n_u=%i, n_a=%i\n",n_klocs,n_u,n_a);
  printf("aInd=%i, %i, %i\n",aInd[0],aInd[1],aInd[2]); 
  printf("AK\n");
  */
  
  for(i=0; i<n_a; i++){
    for(j=0; j<n_u; j++){
      gsl_matrix_set(AK,cntr,j,gsl_matrix_get(K_sigma,aInd[i],j));
      /* printf("%e ",gsl_matrix_get(K_sigma,aInd[i],j));  */
    }
    cntr = cntr + 1;
     /* printf("\n"); */
  } 
      
  /*
  for(i=0; i < n_klocs; i++){ 
    if(i == aInd[cntr_a]){
      for(j=0; j < n_u; j++){ 
        gsl_matrix_set(AK,cntr,j,gsl_matrix_get(K_sigma,i,j));
        printf("%e ",gsl_matrix_get(K_sigma,i,j));
      } 
      printf("\n");
      cntr = cntr + 1;
      cntr_a = cntr_a + 1;
      if(cntr_a == n_a){
        printf("finished\n");
        return;
      }
    }
  }
  printf("\n");
  */
}

/* Intended to multipy Omega_ij * A_ij * K^sigma where A = A_ij * K^sigma
   and assuming Omega_ij is a diagonal precision matrix represented 
   in the values of *v */
void getOmegaAK(gsl_matrix *A, gsl_vector *v, gsl_matrix *B){
  int i,j;
  /* printf("AKOMEGA\n"); */
  for(i=0; i < A->size1; i++){
    for(j=0; j < A->size2; j++){
      gsl_matrix_set(B,i,j,gsl_matrix_get(A,i,j)*gsl_vector_get(v,i));
/*      printf("%e ",gsl_matrix_get(A,i,j)*gsl_vector_get(v,i)); */
    }
/*    printf("\n");*/
  }
  /*printf("\n");*/
}

void getDij(gsl_vector *Dij, gsl_matrix *AK, gsl_vector *prec, gsl_vector *P){
  int i,j;
  double val;
  /*
  printf("\n precs=%e,%e,%e\n",gsl_vector_get(prec,0),gsl_vector_get(prec,1),gsl_vector_get(prec,2));
  printf("ps=%e,%e,%e\n",gsl_vector_get(P,0),gsl_vector_get(P,1),gsl_vector_get(P,2));
  */
  
  for(i=0; i < AK->size2; i++){
    val = 0.0;
    /* printf("i=%i\n",i); */
    for(j=0; j < AK->size1; j++){
      val = val + gsl_vector_get(prec,j) * gsl_vector_get(P,j) * gsl_matrix_get(AK,j,i);
/*       printf("%e ",gsl_vector_get(prec,j) * gsl_vector_get(P,j) * gsl_matrix_get(AK,j,i)); */
    }
/*    printf("\n"); */
    gsl_vector_set(Dij,i,val);
  }
}


double lpost(int N, int ni_max, gsl_vector *v, double tau2_u, double tau2_v, double delta,
            struct data DAT[N][ni_max], int n_u, int n_v, int n_klocs,
            gsl_matrix *K_sigma, gsl_matrix *K_delta,
            double *klocs, gsl_matrix *W, int *NIs, double alpha_u,
            double alpha_v, double beta_u, double beta_v, int draw_u,
            gsl_matrix * U, const gsl_rng * r,int print_warn, int *dpc_status,
            double theC)
{

  double finalterm = 0.0;
  double fterm = 0;
  int i,j,status;
  gsl_matrix * WW = gsl_matrix_alloc(n_u,n_u); /* Make a copy of W */
  if(!WW){
    ERR_MEM_LP
  }
  gsl_matrix_memcpy(WW,W);
  gsl_matrix_scale(WW,1.0/tau2_u); /* NOTE: WW is now W/tau2_u */
  int  status_master = 0;
  #ifdef _OPENMP
    #pragma omp parallel for shared(DAT,K_sigma,r,status_master,U,WW) private(fterm,i,j,status) reduction(+:finalterm) 
  #endif
  for(i=0; i < N; i++){
    gsl_matrix * B = gsl_matrix_alloc(n_u,n_u);
    gsl_vector * Dij = gsl_vector_alloc(n_u);
    gsl_vector * Di = gsl_vector_calloc(n_u);
    gsl_matrix * Ci = gsl_matrix_calloc(n_u,n_u);
    if(!(B && Dij && Di && Ci)){
      status_master = max_int(status_master,2); 
    }
    for(j=0; j < NIs[i]; j++){ 
      /* Allocate Memory */
      gsl_matrix * AK_ij = gsl_matrix_alloc(DAT[i][j].n,n_u);
      gsl_matrix * OmegaAK = gsl_matrix_alloc(DAT[i][j].n,n_u);
      if(!(AK_ij && OmegaAK)){
        status_master = max_int(status_master,2);
      }
      getAK(AK_ij, K_sigma, DAT[i][j].aInd);
      getOmegaAK(AK_ij,DAT[i][j].prec,OmegaAK);
      status = gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,AK_ij,OmegaAK,0.0,B);
      if(status){
        if(print_warn){
          DPC_WARNING(status,gsl_blas_dgemm,lpost());
        }
        status_master = max_int(status_master,2);
      }
      /* Sum the B_ij's into C_i */
      gsl_matrix_add(Ci,B);
      /* Get D_ij */
      getDij(Dij,AK_ij,DAT[i][j].prec,DAT[i][j].p);       
/*
for(kk=0; kk < n_u; kk++){
  printf("%e ",gsl_vector_get(Dij,kk));
}
printf("\n");
*/
      gsl_vector_add(Di,Dij); 
      gsl_matrix_free(AK_ij);
      gsl_matrix_free(OmegaAK);
    }

    /* Add W/tau2_u to Ci; NOTE WW was transformed to W/tau2_u above */
    gsl_matrix_add(Ci,WW);
/*
    for(j=0; j < n_u; j++){
      for(kk=0; kk<n_u; kk++){
        printf("%e ",gsl_matrix_get(Ci,j,kk));
      }
      printf("\n");
    }
    printf("\n");
*/
/*
    for(j=0; j < n_u; j++){
      printf("%e ",gsl_vector_get(Di,j));
    }
*/
    if(draw_u > 0){ /* Draw Us and put into U matrix */
      /* Draw independent random variables into vector, z*/
      gsl_vector * z = gsl_vector_alloc(n_u);
      if(!z){
        status_master = max_int(status_master,2); 
      }
      for(j=0; j < n_u; j++){
        gsl_vector_set(z,j, gsl_ran_ugaussian(r));
      }
      /* Copy the lower diagonal of Ci matrix */

      status = gsl_linalg_cholesky_decomp(Ci); // Ci is now the lower cholesky decomposition
      if(status){
        if(print_warn){
          DPC_WARNING(status,gsl_linalg_cholesky_decomp,lpost());
        }
      }
      if(!status){
        status = gsl_blas_dtrsv(CblasLower,CblasTrans,CblasNonUnit,Ci,z);
        // z is now (C_i^{-1})^T * z)
        if(status){
          if(print_warn){
            DPC_WARNING(status,gsl_blas_dtrsv,lpost());
          }
        }
      }
      if(!status){
        status = gsl_linalg_cholesky_svx(Ci,Di);
        // Di is now Ci^{-1} Di 
        if(status){
          if(print_warn){
            DPC_WARNING(status,gsl_linalg_cholesky_svx,lpost());
          }
        }
      }
      /* finish assigning u */
      if(!status){
        for(j=0; j<n_u; j++){
          gsl_matrix_set(U,i,j,gsl_vector_get(Di,j) + gsl_vector_get(z,j));
        }
      }else{
        for(j=0; j<n_u; j++){
          gsl_matrix_set(U,i,j,-999.0);
        }
      }
      gsl_vector_free(z);
    }else{
      // gsl_matrix_memcpy(Ci2,Ci); // just in case we run into problems
      fterm = get_finalterm(Ci,Di,&status,print_warn);
/*
      if(status == GSL_EDOM){
        gsl_matrix_memcpy(Ci,Ci2);
        // Add noise to diagonal
        for(j=0; j<n_u; j++){
          gsl_matrix_set(Ci,j,j,gsl_matrix_get(Ci,j,j) + \
                                gsl_ran_gaussian_ziggurat(r,10e-9));
        }
        //try again
        fterm = get_finalterm(Ci,Di,&status);
      }
      if(status == GSL_EDOM){
        printf("WARNING: Matrix inversion still failed after adding noise.");
        return(-INFINITY); // 
      }
*/
      if(status){
        status_master = 1;
      }
      finalterm = finalterm + fterm;
    }
      /* WARNING: get_finalterm() modifies Ci to the LU Decomposition */    
    gsl_matrix_free(B);
    gsl_matrix_free(Ci);
    gsl_vector_free(Di);
    gsl_vector_free(Dij);
  }
  if(status_master){
    *dpc_status = status_master;
    return(-INFINITY);
  }
  double vtv;
  gsl_blas_ddot(v,v,&vtv);
  double lpostout = -(.5*N*n_u + 1 + alpha_u) * log(tau2_u) -
                  (.5*n_v + 1 + alpha_v) * log(tau2_v) -
                  vtv/(2.0*tau2_v) -
                  beta_u/tau2_u -
                  beta_v/tau2_v +
                  finalterm + 
                  theC;
/*
  printf("%f\n%f\n%f\n%f\n%f\n%f\n",-(.5*N*n_u + 1 + alpha_u) * log(tau2_u),\
         -(.5*n_v + 1 + alpha_v) * log(tau2_v),\
         -vtv/(2.0*tau2_v),\
         -beta_u/tau2_u,\
         -beta_v/tau2_v,finalterm);
*/
  gsl_matrix_free(WW);
  return(lpostout);
} 

void get_Ksigma(gsl_matrix *K_sigma,int n_klocs, int n_u, double * locs,
                double *k_u, gsl_vector * sigmas){
  int i,j;
  double sigma_t,val, sum;
/*  printf("Ksigma\n"); */
  for(i=0; i < n_klocs; i++){
    sigma_t = gsl_vector_get(sigmas,i);
    sum = 0.0;
    for(j=0; j < n_u; j++){
      val = 1.0/(sqrt(2.0*M_PI)*sigma_t) * 
            exp( - (1.0/(2.0*sigma_t*sigma_t)) * 
                   (locs[i] - k_u[j])*(locs[i] - k_u[j]));
      gsl_matrix_set(K_sigma,i,j,val);
      sum += val;
      /*
      printf("%e ",gsl_matrix_get(K_sigma,i,j));
      if(j == n_u -1){ printf("\n"); } 
      */
    }
    // Normalize row i to ensure the row sum is 1
    for(j=0; j<n_u; j++){
      val = gsl_matrix_get(K_sigma,i,j);
      val /= sum;
      gsl_matrix_set(K_sigma,i,j,val);
    }
  }
}

void get_Kdelta(gsl_matrix *K_delta, int n_klocs, int n_v, double *klocs, double *k_v,
                double delta){
  int i,j;
  double val, sum;
  for(i=0; i<n_klocs; i++){
    sum = 0.0;
    for(j=0; j<n_v; j++){
      val = 1.0/(sqrt(2.0*M_PI)*delta) * 
            exp(-(klocs[i] - k_v[j])*(klocs[i] - k_v[j])/(2.0*delta*delta));
      gsl_matrix_set(K_delta,i,j,val);
      sum += val;
      /* print matrix   */
 /*     printf("%3f ",val);
      if(j == n_v -1){ printf("\n"); }  */
    }
    // Normalize row i to ensure the row sum is 1
    for(j=0; j<n_v; j++){
      val = gsl_matrix_get(K_delta,i,j);
      val /= sum;
      gsl_matrix_set(K_delta,i,j,val);
    }
  }
}

/*
int sum_vector_int(int A[], int n)
{
  int i, sum = 0;
  for(i=0; i < n; i++){
    sum += A[i];
  }
  return sum;
}
*/

/* 
This function modifies the input argument S to be a subset
   of the matrix A.  S = A[ I[0]:I[1] ][ J[0]:j[1] ]
*/
/* 
void submat(int n, int p, int I[2], int J[2], int A[n][p], int S[][J[1]-J[0]+1])
{
  int nout = I[1] - I[0];
  int pout = J[1] - J[0];
  int i,j;
  for(i=I[0]; i < I[1] + 1; i++){
    for(j=J[0]; j < J[1] + 1; j++){
      S[i - I[0]][j - J[0]] = A[i][j];
    }
  }
}
*/

/*
void createindex(int A[]){ */
/*  int * A = realloc(A,sizeof(int)*4); */
/*
  A[0] = 1;
  A[1] = 2;
  A[2] = 3;
  A[3] = 4;
}
*/





