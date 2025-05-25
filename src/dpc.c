//#include "dpc_macros.h"
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include "dpc.h"
#include <string.h>

void dpc(double *k,double *p,double *prec,int *I,int *J,
                int *As, double *klocs, int *NIs, int *N,
                int *ntalldat, double *v_double, int *n_klocs,
                int *n_u, int *n_v, double *k_u, double *k_v,
                double *tau2_u, double *tau2_v, double *delta,
                double *alpha_u, double *alpha_v, double *beta_u,
                double *beta_v,int *nmcmc, int *burn, double *sprop_v,
                double *sprop_tau2_u, double *sprop_delta,
                int *n_print,char **out_dir, int *seed, int *draw_u,
                int *num_threads, double *delta_max,int *dpc_status,
                int *print_warn, int *thin)
{
  gsl_set_error_handler_off(); // Manually check all gsl routines for errors

/*
  if(*print_warn){
    #undef dpc_WARNING
    #define dpc_WARNING(STATUS,F1,F2) \
      printf("WARNING: %s ( in %s ) exited with status %d.\n",\
            #F1,#F2,STATUS);
    #undef dpc_SING_WARNING
    #define dpc_SING_WARNING printf("WARNING: matrix is singular in get_finalterm().\n");
  }else{
    #undef dpc_WARNING
    #define dpc_WARNING(STATUS,F1,F2)
    #undef dpc_SING_WARNING
    #define dpc_SING_WARNING
  }
*/


  #ifdef _OPENMP
    #pragma omp parallel
    {
      if(*num_threads == 0){
        *num_threads = omp_get_num_threads();
      }else if(*num_threads > omp_get_num_threads()){
        *num_threads = omp_get_num_threads();
        #pragma omp master
        printf("NOTE: Number of requested threads > available threads.\n");
      }
      omp_set_num_threads(*num_threads);
      #pragma omp master
      printf("NOTE: %d/%d threads used.\n",*num_threads,\
             omp_get_num_threads());
    }
  #else
    printf("OpenMP not available.\n");
  #endif

  /*
      k: list of all locations
      p: list of all response variables
      I: Index on which smooth process k,p compe from
      J: Index on which replication from the ith process
      NIs: The number of reps for each i
      N: The number of smooth functions, (indexed by i) we are estimating
      max_ni: largest value of ni
  */
  /* Create an index indicating where each distinct dataset is located */
  /* Find the largest number of reps for a mean function */
  int i,j,kk,ok;

  int ni_max = 0;
  for(i=0; i < *N; i++){
    if(NIs[i] > ni_max) ni_max = NIs[i];
  }

  struct data DAT[*N][ni_max];
  /* format data in DAT */
  data_creator(As,I,J,k,*N,NIs,ni_max,*ntalldat,p,prec,DAT,dpc_status);
  if(*dpc_status) return;

  // Obtain a portion of the constant of the log posterior
  double theC;
  theC = get_C(*N,NIs,ni_max,DAT);
  //  printf("theC=%.*e\n",MYE,theC);

  /* INITIAL VALUES */
  /* Create W (BROWNIAN MOTION) */
  gsl_matrix * W = gsl_matrix_alloc(*n_u,*n_u);
  if(!W){
    ERR_MEM
  }
  brownian_precision(W);

  /* Put V in gsl format */
  gsl_vector * v = gsl_vector_alloc(*n_v);
  if(!v){
    ERR_MEM
  }
  for(i=0; i<*n_v; i++){
    gsl_vector_set(v,i,v_double[i]);
  }

  /* Get Starting Values */
  gsl_matrix * K_sigma = gsl_matrix_alloc(*n_klocs,*n_u);
  gsl_matrix * K_delta = gsl_matrix_alloc(*n_klocs,*n_v);
  gsl_matrix * U = gsl_matrix_alloc(*N,*n_u);
  if(!(K_sigma && K_delta && U)){
    ERR_MEM
  }
  get_Kdelta(K_delta,*n_klocs,*n_v,klocs,k_v,*delta);
  gsl_vector * sigmas = gsl_vector_alloc(*n_klocs);
  if(!sigmas){
    ERR_MEM
  }
  gsl_blas_dgemv(CblasNoTrans,1.0,K_delta,v,0.0,sigmas);

  // Check to make sure all sigmas are positive...
  *dpc_status = 0;
  ok=1;
  for(i=0; i<*n_klocs; i++){
    if(gsl_vector_get(sigmas,i) <= 0){
      ok=0;
      break;
    }
  }
  if(!ok){
    printf("ERROR: Starting delta and v values resulted in negative sigmas. Please use different starting values.\n");
    *dpc_status = 1;
    return;
  }


  /*
  for(i=0; i<*n_u; i++){
    printf("%f\n",gsl_vector_get(sigmas,i));
  }
  */
  get_Ksigma(K_sigma, *n_klocs, *n_u, klocs, k_u, sigmas);


  /* INITIALIZE RANDOM SEED */

  const gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
  if(!r){
    ERR_MEM
  }
  unsigned long int theseed;

  if(*seed == 0){
    theseed = time(NULL);
  }else{
    theseed = *seed;
  }

  gsl_rng_set(r,theseed);

  /* Write a log likelihood function */
  double lf;
  lf = lpost(*N,ni_max,v,*tau2_u, *tau2_v,*delta, DAT, *n_u, *n_v,
        *n_klocs, K_sigma, K_delta, klocs, W, NIs, *alpha_u,
        *alpha_v, *beta_u, *beta_v,0,U,r,*print_warn,dpc_status,theC);
  if(isinf(lf)){
    printf("ERROR: First posterior evaluation was non-finite. Please use different starting values and/or grid.");
    *dpc_status = 1;
    return;
  }

  // printf("Initial lf=%15f\n",lf);  


  /* Initialize values */
  gsl_vector * vstar = gsl_vector_alloc(*n_v);
  gsl_vector * sigmas_star = gsl_vector_calloc(*n_klocs);
  gsl_vector * cntr_v = gsl_vector_calloc(*n_v);
  gsl_matrix * K_sigma_star = gsl_matrix_alloc(*n_klocs,*n_u);
  gsl_matrix * K_delta_star = gsl_matrix_alloc(*n_klocs, *n_v);
  gsl_matrix * UU = gsl_matrix_alloc(*nmcmc * *N, *n_u + 2);
  if(!(vstar && sigmas_star && cntr_v && K_sigma_star && K_delta_star && UU)){
    ERR_MEM
  }
  double lfstar, lr;
  double delta_star, tau2_u_star;
  int cntr_delta = 0;
  int cntr_tau2_u = 0;
  int cntr_out = 0;
  int status=0;
  /* Create vectors for storage */
  double * TAU2_U = (double *)malloc(sizeof(double)* *nmcmc);
  double * TAU2_V = (double *)malloc(sizeof(double)* *nmcmc);
  double * DELTA = (double *)malloc(sizeof(double)* *nmcmc);
  double * V = (double *)malloc(sizeof(double)* *nmcmc * *n_v);
  double * LLIKE = (double *)malloc(sizeof(double) * *nmcmc);

  if(!(TAU2_U && TAU2_V && DELTA && V && LLIKE)){
    ERR_MEM
  }
  /* START MCMC */
  time_t end, begin;
  time(&begin);
  double time_spent, time_remaining;
  int iter, burn_num;
  gsl_vector_memcpy(vstar,v);
  printf("Starting MCMC...\n");
  for(i=0; i<(*nmcmc * *thin + *burn); i++){
    /* Draw vs */
    for(j=0; j<*n_v; j++){
      /* Draw a proposal */
      gsl_vector_set(vstar,j,gsl_vector_get(v,j) + gsl_ran_gaussian(r,sprop_v[j]));
      /* Update sigmas based on v */
      status = gsl_blas_dgemv(CblasNoTrans,1.0,K_delta,vstar,0.0,sigmas_star);
      if(status){
        if(*print_warn){
          DPC_WARNING(status,gsl_blas_dgemv,dpc());
        }
        ok=0;
      }else{
        ok=1;
      }
      /* make sure sigmas are all positive */
      if(ok){
        for(kk=0; kk < *n_klocs; kk++){
          if(gsl_vector_get(sigmas_star,kk) <= 0){
            ok = 0;
            break;
          }
        }
      }
      if(ok){ // All sigmas are positive
        /* Update Sigma_star */
        get_Ksigma(K_sigma_star, *n_klocs, *n_u, klocs, k_u,sigmas_star);
        /* Calculate lpost */
        lfstar = lpost(*N,ni_max,vstar,*tau2_u, *tau2_v,*delta, DAT, *n_u, *n_v,
           *n_klocs, K_sigma_star, K_delta, klocs, W, NIs, *alpha_u,
           *alpha_v, *beta_u, *beta_v,0,U,r,*print_warn,dpc_status,theC);
        if(*dpc_status) return;
//printf("v%d: %f  ...",j,exp(lfstar-lf));
        lr = lfstar - lf;
        /* Accept or Reject (MH step) */
        if(lr > log(gsl_ran_flat(r,0,1))){
//printf("accepted");
          gsl_vector_set(v,j,gsl_vector_get(vstar,j));
          lf = lfstar;
          gsl_matrix_memcpy(K_sigma,K_sigma_star);
          gsl_vector_set(cntr_v,j,gsl_vector_get(cntr_v,j) + 1);
        }else{
          gsl_vector_set(vstar,j,gsl_vector_get(v,j));
        }
      }else{ // If not all sigmas are positive, reset vstar
        gsl_vector_set(vstar,j,gsl_vector_get(v,j));
      }
//printf("\n");
    }
    /* Sample tau_u */
    // the proposal standard deviation is the smaller of 
    //    the sd given and the distance to 0 of tau2_u
//    tau2_u_star = rtruncnorm(r,0,INFINITY,*tau2_u,min_double(*sprop_tau2_u,*tau2_u));
    tau2_u_star = rtruncnorm(r,0,INFINITY,*tau2_u,*sprop_tau2_u);
    lfstar = lpost(*N,ni_max,v,tau2_u_star, *tau2_v,*delta, DAT, *n_u, *n_v,
         *n_klocs, K_sigma, K_delta, klocs, W, NIs, *alpha_u,
         *alpha_v, *beta_u, *beta_v,0,U,r,*print_warn,dpc_status,theC);
    if(*dpc_status) return;
    lr = lfstar - lf + log(dtruncnorm(*tau2_u,0,INFINITY,tau2_u_star,*sprop_tau2_u)) -
                       log(dtruncnorm(tau2_u_star,0,INFINITY,*tau2_u,*sprop_tau2_u));
//printf("tau2_u pr=%f, r=%f...",exp(lfstar-lf),exp(lr));
    if(lr > log(gsl_ran_flat(r,0,1))){
//printf("accepted");
      *tau2_u = tau2_u_star;
      lf = lfstar;
      cntr_tau2_u = cntr_tau2_u + 1;
    }
//printf("\n");
    /* Sample tau2_v */
    *tau2_v = 1/gsl_ran_gamma(r,*alpha_v + *n_v/2.0,1.0/ (.5*sum2(v) + *beta_v));
    lf = lpost(*N,ni_max,v,*tau2_u, *tau2_v,*delta, DAT, *n_u, *n_v,
         *n_klocs, K_sigma, K_delta, klocs, W, NIs, *alpha_u,
         *alpha_v, *beta_u, *beta_v,0,U,r,*print_warn,dpc_status,theC);
    if(*dpc_status) return;
    /* Sample Delta */
    // The proposal sd is the smaller of *sprop_delta and *delta's distance to a boundary
/*    
    delta_star = rtruncnorm(r,0,*delta_max,*delta,\
       min_double(*sprop_delta,min_double(*delta,*delta_max - *delta)));
*/
    delta_star = rtruncnorm(r,0,*delta_max,*delta,*sprop_delta);
// printf("delta=%f, delta_star=%f, sigma=%e\n",*delta,delta_star,min_double(*sprop_delta,min_double(*delta,*delta_max - *delta)));
    get_Kdelta(K_delta_star,*n_klocs,*n_v,klocs,k_v,delta_star);
    /* Make sure delta yields positive sigmas */
    status = gsl_blas_dgemv(CblasNoTrans,1.0,K_delta_star,v,0.0,sigmas_star);
    if(status){
      if(*print_warn){
        DPC_WARNING(status,gsl_blas_dgemv,dpc());
      }
      ok = 0;
    }else{
      ok = 1;
    }
    /* make sure sigmas are all positive */
    if(ok){
      for(kk=0; kk < *n_klocs; kk++){
        if(gsl_vector_get(sigmas_star,kk) <= 0){
          ok = 0;
          break;
        }
      }
    }
    if(ok){
      get_Ksigma(K_sigma_star, *n_klocs, *n_u, klocs, k_u, sigmas_star);
      lfstar = lpost(*N,ni_max,v,*tau2_u, *tau2_v,delta_star, DAT, *n_u, *n_v,
           *n_klocs, K_sigma_star, K_delta_star, klocs, W, NIs, *alpha_u,
           *alpha_v, *beta_u, *beta_v,0,U,r,*print_warn,dpc_status,theC);
      if(*dpc_status) return;
      lr = lfstar - lf + log(dtruncnorm(*delta,0,*delta_max,delta_star,*sprop_delta)) -
                       log(dtruncnorm(delta_star,0,*delta_max,*delta,*sprop_delta));

      if(lr > log(gsl_ran_flat(r,0,1))){
//printf("accepted");
        *delta = delta_star;
        lf = lfstar;
        cntr_delta = cntr_delta + 1;
        gsl_matrix_memcpy(K_delta,K_delta_star);
        gsl_matrix_memcpy(K_sigma,K_sigma_star);
      }
//printf("\n");
    }// If not ok, it is automatically rejected, no need to do anything else...
    /* Draw U */
    if(*draw_u > 0){
      lpost(*N,ni_max,v,*tau2_u, *tau2_v,*delta, DAT, *n_u, *n_v,
           *n_klocs, K_sigma, K_delta, klocs, W, NIs, *alpha_u,
           *alpha_v, *beta_u, *beta_v,1,U,r,*print_warn,dpc_status,theC);
      if(*dpc_status) return;
    }

    /* Print Progress */
    time(&end);
    time_spent = difftime(end,begin)/60.0;
    iter = i + 1 - *burn;
    if(iter < 0){ iter = 0;}
    burn_num = i + 1;
    if(burn_num > *burn){ burn_num = *burn;}
    if(((i+1) % *n_print) == 0){
      printf("iter=%i, burn=%i, log-like=%.2f, minutes=%.2f\n",
             iter,burn_num,lf,time_spent);
      printf("    acceptances: tau2_u = %.2f, delta = %.2f\n    ",
          cntr_tau2_u/(1.0 * i + 1.0),cntr_delta/(1.0 * i + 1.0));
      for(j=0; j<*n_v; j++){
        printf("v%d = %.2f",j+1,gsl_vector_get(cntr_v,j)/(1.0*i + 1.0));
        if(j < (*n_v - 1)){
          printf(", ");
        }else{
          printf("\n");
        }
      }

      time_remaining = time_spent/(i*1.0 + 1.0) * ( (*nmcmc * *thin + *burn) - (i + 1) );
      printf("Minutes Remaining: %.2f\n\n",time_remaining);

    }
    /* Save samples */
    if( ((i + 1) > *burn) & (((i + 1) % *thin) == 0)){
      // cntr_out = i - *burn;	  
      TAU2_U[cntr_out] = *tau2_u;
      TAU2_V[cntr_out] = *tau2_v;
      DELTA[cntr_out] = *delta;
      for(j=0; j<*n_v; j++){
        // *(V + cntr_out * *n_v + j) = gsl_vector_get(v,j);
        V[cntr_out* *n_v + j] = gsl_vector_get(v,j);
      }
      LLIKE[cntr_out] = lf;
      /* Store U's */
      if(*draw_u){
        for(j=0; j < *N; j++){
          gsl_matrix_set(UU,*N * cntr_out + j,0,cntr_out + 1);
          gsl_matrix_set(UU,*N * cntr_out + j,1, j + 1);
          for(kk=0; kk<*n_u; kk++){
            gsl_matrix_set(UU,*N * cntr_out + j,kk + 2, gsl_matrix_get(U,j,kk));
          }
        }
      }
	  cntr_out = cntr_out + 1;
    }
  }
//  end = clock();
//  clock_gettime(CLOCK_MONOTONIC,&end);
  time(&end);
  time_spent = difftime(end,begin)/60.0;
//  time_spent =(double)(end.tv_sec - begin.tv_sec)/(60.0);

  printf("MCMC Complete!\n\n");

  /* Write output to file */
  char * out_dir2 = malloc(1000);
  if(!out_dir2){
    ERR_MEM
  }

  /* Write Info File */
  sprintf(out_dir2,"%s%s",*out_dir,"mcmc_info.txt");
  FILE *ifile = fopen(out_dir2,"w");
  fprintf(ifile,"Total minutes     = %f\n",time_spent);
  fprintf(ifile,"MCMC collected    = %d\nBurn in           = %d\n",*nmcmc,*burn);
  fprintf(ifile,"Thinning          = %d\n",*thin);
  fprintf(ifile,"tau2_u acceptance = %f\n",cntr_tau2_u/(1.0* *nmcmc * *thin + *burn));
  fprintf(ifile,"delta  acceptance = %f\n",cntr_delta/(1.0 * *nmcmc * *thin + *burn));
  for(i=0; i < *n_v; i++){
    fprintf(ifile,"v%02d  acceptance = %f\n",i+1,gsl_vector_get(cntr_v,i)/(1.0**nmcmc * *thin + *burn));
  }
  fclose(ifile);



  sprintf(out_dir2,"%s%s",*out_dir,"output.csv");
  FILE *f = fopen(out_dir2,"w");
  sprintf(out_dir2,"%s%s",*out_dir,"v_output.csv");
  FILE *f2 = fopen(out_dir2,"w");
  if(f == NULL || f2 == NULL){
    printf("Error opening file! Results cannot be saved. :(\n");
  }
  fprintf(f,"%s\n","tau2_u,tau2_v,delta");
  for(i=0; i<*nmcmc; i++){
    fprintf(f,"%.*e, %.*e, %.*e\n",MYE,TAU2_U[i],MYE,TAU2_V[i],MYE,DELTA[i]);
    for(j=0; j<(*n_v - 1); j++){
      fprintf(f2,"%.*e,",MYE,*(V + i * *n_v + j));
    }
    fprintf(f2,"%.*e\n",MYE,*(V + i * *n_v + j));
  }
  fclose(f);
  fclose(f2);
  /* Write log-like, k_u, k_v, and the locations to file */
  sprintf(out_dir2,"%s%s",*out_dir,"k_u.txt");
  FILE * ufile = fopen(out_dir2,"w");
  sprintf(out_dir2,"%s%s",*out_dir,"k_v.txt");
  FILE * vfile = fopen(out_dir2,"w");
  sprintf(out_dir2,"%s%s",*out_dir,"k_locs.txt");
  FILE * lfile = fopen(out_dir2,"w");
  sprintf(out_dir2,"%s%s",*out_dir,"log_post.txt");
  FILE * LLfile = fopen(out_dir2,"w");
  if(ufile == NULL || vfile == NULL || lfile == NULL || LLfile == NULL){
    printf("Error opening file! Results cannot be saved (second group). :(\n");
  }
  /* k_u */
  for(i=0; i<(*n_u-1); i++){
    fprintf(ufile,"%.*e\n",MYE,k_u[i]);
  }
  fprintf(ufile,"%.*e",MYE,k_u[i]);
  /* k_v */
  for(i=0; i<(*n_v - 1); i++){
    fprintf(vfile,"%.*e\n",MYE,k_v[i]);
  }
  fprintf(vfile,"%.*e",MYE,k_v[i]);
  /* k_locs */
  for(i=0; i<(*n_klocs -1); i++){
    fprintf(lfile,"%.*e\n",MYE,klocs[i]);
  }
  fprintf(lfile,"%.*e",MYE,klocs[i]);
  for(i=0; i<(*nmcmc -1); i++){
    fprintf(LLfile,"%.*e\n",MYE,LLIKE[i]);
  }
  fprintf(LLfile,"%.*e",MYE,LLIKE[i]);
  fclose(ufile);
  fclose(vfile);
  fclose(lfile);
  fclose(LLfile);
  /* Write U's to file */
/*
  if(*draw_u > 0){
    sprintf(out_dir2,"%s%s",*out_dir,"Us.csv");
    FILE * Ufile = fopen(out_dir2,"w");
    if(Ufile == NULL){
      printf("Error opening U file.  Results cannot be saved. :(\n");
    }
    for(i=0; i < UU->size1; i++){
      for(j=0; j < (UU->size2 - 1); j++){
        fprintf(Ufile,"%f,",gsl_matrix_get(UU,i,j));
      }
      fprintf(Ufile,"%f\n",gsl_matrix_get(UU,i,j));
    }
    fclose(Ufile);
  }
*/
  /* Write U's to file (alternate -- multiple files) */
//  FILE * Ufile;
  if(*draw_u > 0){
    for(i=0; i < *N; i++){
      sprintf(out_dir2,"%s%s%04d%s",*out_dir,"us/u",i+1,".csv");
      FILE * Ufile = fopen(out_dir2,"w");
      if(Ufile == NULL){
        printf("Error opening file to save U%d.  Results cannot be saved. :(\n",i);
      }
      for(j=0; j < *nmcmc; j++){
        for(kk=2; kk < (*n_u+1); kk++){
          fprintf(Ufile,"%.*e,",MYE,gsl_matrix_get(UU,*N * j + i, kk));
        }
        fprintf(Ufile,"%.*e\n",MYE,gsl_matrix_get(UU,*N * j + i, kk));
      }
      fclose(Ufile);
    }
  }
  gsl_vector_free(cntr_v);
  gsl_vector_free(sigmas_star);
  gsl_vector_free(vstar);
  free(out_dir2);
  free(TAU2_U);
  free(TAU2_V);
  free(DELTA);
  free(V);
  free(LLIKE);
  gsl_matrix_free(K_sigma_star);
  gsl_matrix_free(K_delta_star);
  /*
  printf("\nlf=%e\n",lf);
  */

  /*
  int * A = malloc(4*sizeof(int));
  createindex(A);
  printf("%3i %3i %3i %3i\n",A[0],A[1],A[2],A[3]);
  free(A);
  */


  /* Create an index on k, p */
  /*
  int num_rows = sizeof(p) / sizeof(p[0]);
  printf("The number of rows is: %6i \n",num_rows);
  printf("The first value is: %6f \n",p[1]);
  printf("The number of rows is: %6i \n",sum_vector_int(NIs,*N));
  int idxptr[*N][ni_max];
  for(i=0; i<*N; i++){
    for(j=0; j<NIs[i]; j++){
      for(kk=0; kk<*ntalldat; kk++){
        if(I[kk] == i && J[kk] == j){
        }
      }
    }
  }
  double *b = k + 130;
  printf("%3f %3f \n", b[1],b[2]);
  */
  gsl_matrix_free(K_sigma);
  gsl_matrix_free(K_delta);
  gsl_matrix_free(W);
  gsl_matrix_free(U);
  gsl_matrix_free(UU);
  gsl_vector_free(sigmas);
  gsl_vector_free(v);
  for(i=0; i < *N; i++){
    for(j=0; j < NIs[i]; j++){
      free(DAT[i][j].aInd);
    }
  }
}




void udraw(double * alpha_u, double * alpha_v, int *As,
           double * beta_u, double * beta_v,
           double * delta, int * I, int * J, double * k,
           double * k_u, double * k_v, double * locs,
           int *N, int * n_klocs, int * n_u, int * n_v,
           int * NIs, int * ni_max, int * nmcmc,
           int *ntalldat, double *p, double *prec, int * seed,
           double * tau2_u, double * tau2_v, double * UU,
           double * V, int *print_warn,int *dpc_status){
/*
  if(*print_warn){
    #undef DPC_WARNING
    #define DPC_WARNING(STATUS,F1,F2) \
      printf("WARNING: %s ( in %s ) exited with status %d.\n",\
            #F1,#F2,STATUS);
  }else{
    #undef DPC_WARNING
    #define DPC_WARNING(STATUS,F1,F2)
  }
*/
  const gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
  if(!r){
    ERR_MEM
  }
  unsigned long int theseed;
  if(*seed == 0){
    theseed = time(NULL);
  }else{
    theseed = *seed;
  }
  gsl_rng_set(r,theseed);

  /* Format Data */
  struct data DAT[*N][*ni_max];
  data_creator(As,I,J,k,*N,NIs,*ni_max,*ntalldat,p,prec,DAT,dpc_status);
  if(*dpc_status) return;


  gsl_matrix * K_delta = gsl_matrix_alloc(*n_klocs,*n_v);
  gsl_matrix * K_sigma = gsl_matrix_alloc(*n_klocs,*n_u);
  gsl_matrix * U = gsl_matrix_alloc(*N,*n_u);
  gsl_matrix * W = gsl_matrix_alloc(*n_u,*n_u);
  gsl_vector * sigmas = gsl_vector_alloc(*n_klocs);
  gsl_vector * v = gsl_vector_alloc(*n_v);
  if(!(K_delta && K_sigma && U && W && sigmas && v)){
    ERR_MEM
  }

  brownian_precision(W);
  double theC;
  theC = get_C(*N,NIs,*ni_max,DAT);

  int i,j,kk,status;
  int cntr = 0;
  for(i=0; i < *nmcmc; i++){
//printf("i=%d, ",i);
    for(j=0; j<*n_v; j++){
      gsl_vector_set(v,j,V[j* *nmcmc + i]);
    }
    get_Kdelta(K_delta, *n_klocs, *n_v, locs, k_v, delta[i]);
    status = gsl_blas_dgemv(CblasNoTrans,1.0,K_delta,v,0.0,sigmas);
//printf("status=%d",status);
    if(status){
      if(*print_warn){
        DPC_WARNING(status,gsl_blas_dgemv,udraw());
      }
    }else{
      get_Ksigma(K_sigma, *n_klocs, *n_u, locs, k_u, sigmas);

      lpost(*N,*ni_max,v,tau2_u[i],tau2_v[i],delta[i],DAT,*n_u,*n_v,
            *n_klocs, K_sigma, K_delta, locs, W, NIs, *alpha_u,
            *alpha_v, *beta_u, *beta_v,1,U,r,*print_warn,dpc_status,theC);
      if(*dpc_status) return;
    }
    for(j=0; j<*N; j++){
      for(kk=0;kk< (*n_u + 2); kk++){
        if(kk == 0){
          UU[cntr] = i + 1;
        }else if(kk == 1){
          UU[cntr] = j + 1;
        }else{
          if(!status){
            UU[cntr] = gsl_matrix_get(U,j,kk-2);
          }else{
            UU[cntr] = -999.0;
          }
        }
        cntr++;
      }
    }
  }
  gsl_matrix_free(K_delta);
  gsl_matrix_free(K_sigma);
  gsl_vector_free(sigmas);
  gsl_matrix_free(U);
  gsl_vector_free(v);
  gsl_matrix_free(W);
}



/* Function to speed up drawing figures */
void postdraws(double *A, double *delta, int *funcnums, double *k_u, double *k_v,
               double *locs, int *nfunc, int *nfunc_tot, int *ngrid, double *gridvals,
               int *nmcmc,
               int *n_klocs, int *n_u, int *n_v,
               double *U, double *vs,int *print_warn,int *dpc_status){
/*
  if(*print_warn){
    #undef DPC_WARNING
    #define DPC_WARNING(STATUS,F1,F2) \
      printf("WARNING: %s ( in %s ) exited with status %d.\n",\
            #F1,#F2,STATUS);
  }else{
    #undef DPC_WARNING
    #define DPC_WARNING(STATUS,F1,F2)
  }
*/
  gsl_matrix * K_delta = gsl_matrix_alloc(*ngrid,*n_v);
  gsl_matrix * K_sigma = gsl_matrix_alloc(*ngrid,*n_u);
  gsl_vector * sigmas = gsl_vector_alloc(*ngrid);
  gsl_vector * v = gsl_vector_alloc(*n_v);
  if(!(K_delta && K_sigma && sigmas && v)){
    ERR_MEM
  }
  int i,ii,j,k,ok,status;
  double sum;

  // Create Grid
/*
  double gridint;
  double gridvals[*ngrid];
  gridint = (locs[*n_klocs - 1] - locs[0]) / (1.0* *ngrid - 1.0);
  gridvals[0] = locs[0];
  for(i=1; i<*ngrid; i++){
    gridvals[i] = gridvals[i-1] + gridint;
  }
*/

  for(i=0; i<*nmcmc; i++){
    for(j=0; j<*n_v; j++){
      gsl_vector_set(v,j,vs[i + j * *nmcmc]);
    }
    get_Kdelta(K_delta,*ngrid,*n_v,gridvals,k_v,delta[i]);
    status = gsl_blas_dgemv(CblasNoTrans,1.0,K_delta,v,0.0,sigmas);
    if(status){
      if(*print_warn){
        DPC_WARNING(status,gsl_blas_dgemv,postrdraws());
      }
      ok = 0;
    }else{
      ok = 1;
    }
    if(ok){
      get_Ksigma(K_sigma,*ngrid,*n_u,gridvals,k_u,sigmas);
      // Make sure all sigmas are positive
      for(j=0; j<*ngrid; j++){
        if(gsl_vector_get(sigmas,j) < 0){
          ok = 0;
          break;
        }
      }
    }
    if(ok){
      for(j=0; j<*nfunc; j++){
        // Multiply Ksigma and U
        for(ii=0; ii<*ngrid; ii++){
          sum = 0.0;
          for(k=0; k<*n_u; k++){
            sum = sum + gsl_matrix_get(K_sigma,ii,k) * U[i * *nfunc_tot + (funcnums[j]-1) + k * *nmcmc * *nfunc_tot];
          }
          A[j * *nmcmc * *ngrid + i* *ngrid + ii] = sum;
        }
      }
    }else{
      for(j=0; j<*nfunc; j++){
        for(ii=0; ii<*ngrid; ii++){
          A[j * *nmcmc * *ngrid + i* *ngrid + ii] = -999;
        }
      }
    }
  }


  gsl_matrix_free(K_delta);
  gsl_matrix_free(K_sigma);
  gsl_vector_free(sigmas);
  gsl_vector_free(v);
}
/*
void threedarray(int * A){
  int i,j,k;
  int n,p,q;
  int cntr = 1;
  n=2;
  p=3;
  q=4;

  cntr = 1;
  for(k=0; k<q; k++){
    for(j=0; j<p; j++){
      for(i=0; i<n; i++){
        A[k*n*p + j*n + i] = cntr;
        cntr ++;
      }
    }
  }
}
*/




