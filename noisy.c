#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>

gsl_rng * r;

void single_simul_obs(int* y1, int* y2,
		      int i,
		      double tau,
		      double* true_theta,
		      double* prob) {
  prob[0] = true_theta[0] * y1[i-1] * tau;
  //  printf("%f\t", prob[0]);
  prob[1] = prob[0] + true_theta[1] * y1[i-1] * y2[i-1] * tau;
  //  printf("%f\t", prob[1]);
  prob[2] = prob[1] + true_theta[2] * y2[i-1] * tau;
  //  printf("%f\n", prob[2]);
  double u = gsl_rng_uniform(r);
  //printf("\n%f\n\n", u);
  if (u <= prob[0]){
    y1[i] = y1[i-1] + 1;
    y2[i] = y2[i-1];
  } else if (u > prob[0] && u <= prob[1] && y1[i - 1] > 0) {
    y1[i] = y1[i-1] - 1;
    y2[i] = y2[i-1] + 1;
  } else if (u > prob[1] && u <= prob[2] && y2[i - 1] > 0) {
    y1[i] = y1[i-1];
    y2[i] = y2[i-1] - 1;
  } else {
    y1[i] = y1[i-1];
    y2[i] = y2[i-1];
  } 
}
   
void simulation_obs (int* y1, int* y2,
		     double* true_theta,
		     double* prob,
		     int n_obs,
		     double tau) {

  for (int i = 1; i < n_obs; i++)
    single_simul_obs(y1,y2,i,tau,true_theta,prob);
}

double kernel (double true_y1,
	       double true_y2,
	       double sim_y1,
	       double sim_y2,
	       double toler) {
  double res =  gsl_ran_bivariate_gaussian_pdf((true_y1 - sim_y1) / toler,
					       (true_y2 - sim_y2) / toler,
					       100, 100, 0);
  //printf("%f\n",res);
  return res;
}

void prior_theta (double** theta_part,
		  int n_part) {
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < n_part; j++){
      theta_part[i][j] = gsl_rng_uniform(r) * 0.6;
      //printf("%f\n",theta_part[i][j]);
    }
}


void seq_abc (int n_obs,
	      int n_part,
	      int* y1,
	      int* y2,
	      double* true_theta,
	      double** theta_part,
	      double tau,
	      double* output_theta,
	      int noisy) {

  double prob[3];
  double toler = sqrt(tau);
  double x1, x2;
  double w[n_part];
  double prob_theta[n_part + 1];
  int* sim_y1[n_part];
  int* sim_y2[n_part];
  for (int i = 0; i < n_part; i++){
    sim_y1[i] = malloc(n_obs * sizeof(int));
    sim_y2[i] = malloc(n_obs * sizeof(int));
  }

  // First "simulated" obs are true obs
  for (int i = 0; i < n_part; i++){
    sim_y1[i][0] = y1[0];
    sim_y2[i][0] = y2[0];
  }
  
  for (int j = 1; j < n_obs; j++){

    // Step 1
    for (int i = 0; i < n_part; i++){
      single_simul_obs(sim_y1[i],sim_y2[i],j,tau,theta_part[i],prob);
      //if (j == 20) printf("%d\t%d\n", sim_y1[i][j],sim_y2[i][j]);
    }
    
    // Step 2
    prob_theta[0] = 0;
    if (noisy == 0){
      gsl_ran_bivariate_gaussian(r, 0.1, 0.1, 0, &x1, &x2);
      for (int i = 0; i < n_part; i++){
	w[i] = kernel(y1[j], y2[j], sim_y1[i][j] + toler * x1,
		      sim_y2[i][j] + toler * x2, toler);
	prob_theta[i+1] = prob_theta[i] + w[i];
      }
    } else {
      for (int i = 0; i < n_part; i++){
	w[i] = kernel(y1[j], y2[j], sim_y1[i][j], sim_y2[i][j], toler);
	prob_theta[i+1] = prob_theta[i] + w[i];
	//printf("\n%f\n",w[i]);
	//if(j==1) printf("%f\n",prob_theta[i+1]);
      }
    }
    
    // Step 3
    for (int i = 0; i < n_part; i++){
      int flag = -1;
      int l = 0;
      double u = gsl_rng_uniform(r) * prob_theta[n_part];
      while (flag == -1){
	if (u <= prob[l+1] && u > prob[l]) flag = 1;
	l++;
      }
      //printf("\n%f\t%f\t%d\n",prob_theta[n_part],u,l-1);
      for (int k = 0; k < 3; k++)
	theta_part[i][k] = theta_part[l-1][k];
      //printf("%f\t%f\t%f\n",theta_part[i][0],theta_part[i][1],theta_part[i][2]);
    }
  }

  // Output theta
  double u = gsl_rng_uniform(r) * prob_theta[n_part];
  int flag = -1;
  int l = 0;
  while (flag == -1){
    if (u <= prob[l+1] && u > prob[l]) flag = 1;
    l++;
  }
  //printf("\n%f\t%d\n",u,l-1);
  for (int k = 0; k < 3; k++)
    output_theta[k] = theta_part[l-1][k];
  if (noisy == 0)
    printf("\nEstimated parameters ABC with observations:%d\n%f\t%f\t%f\n",n_obs,output_theta[0],output_theta[1],output_theta[2]);
  else
    printf("\nEstimated parameters noisy ABC with observations:%d\n%f\t%f\t%f\n",n_obs,output_theta[0],output_theta[1],output_theta[2]); 
}

int main () {

  r = gsl_rng_alloc(gsl_rng_default);
  double true_theta[3] = {0.5, 0.0025, 0.3};
  double output_theta[3];
  int n_obs = 100;
  double tau = 0.05;
  double prob[3];
  int n_part = 30;
  double h = sqrt(tau);
  int* y1 = malloc(sizeof(int)*n_obs);
  int* y2 = malloc(sizeof(int)*n_obs);

  double *theta_part[n_part];
  for (int i = 0; i < n_part; i++)
    theta_part[i] = malloc(sizeof(double)*3);

  y1[0] = 10;
  y2[0] = 10;  

  printf("\nTrue parameters:\n%f\t%f\t%f\n",true_theta[0],true_theta[1],true_theta[2]);

  // Initialisazion
  simulation_obs(y1, y2, true_theta, prob, n_obs, tau);
  prior_theta(theta_part, n_part);

  //for(int i = 0; i < n_obs; i++)
  //  printf("%d\t%d\n", y1[i], y2[i]);
  
  // Sequential ABC
  seq_abc(n_obs,n_part,y1,y2,true_theta,theta_part,tau,output_theta,1);

  // Sequential noisy-ABC
  seq_abc(n_obs,n_part,y1,y2,true_theta,theta_part,tau,output_theta,0);

  return 0;
}
