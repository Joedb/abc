#include <stdio.h>
#include <gls/gsl_rng.h>
#include <math.h>

const gsl_rng_type * T;
gsl_rng * r;
T = gsl_rng_default;
r = gsl_rng_alloc (T);

void single_simul_obs(int *y1, int *y2,
		      int i;
		      double tau,
		      double *true_theta,
		      double *prob) {
  prob[0] = true_theta[0]*y1[i-1]*tau;
  prob[1] = prob[0] + true_theta[1]*y1[i-1]*y2[i-1]*tau;
  prob[2] = prob[1] + true_theta[2]*y2[i-1]*tau;
  double u = gsl_rng_uniform(r);
  if (u <= prob[0]){
    y1[i] = y1[i-1] + 1;
    y2[i] = y2[i-1];
  } else if (u > prob[0] && u <= prob[1] && y1[i - 1] > 0){
    y1[i] = y1[i-1] - 1;
    y2[i] = y2[i-1] + 1;
  } else if (u > prob[1] && u <= prob[2] && y2[i - 1] > 0){
    y1[i] = y1[i-1];
    y2[i] = y2[i-1] - 1;
  } else {
    y1[i] = y1[i-1];
    y2[i] = y2[i-1];
  } 
}
   
void simulation_obs(int *y1, int *y2,
		    double *true_theta,
		    double *prob,
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
  return gsl_ran_bivariate_gaussian_pdf((true_y1 - sim_y1) / toler,
					(true_y2 - simy_2) / toler,
					1, 1, 0);
}

void prior_theta (double **theta_part,
		  int n_part) {
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < n_part; j++)
      theta_part[i][j] = gls_rng_uniform(r) * 0.6;
}


void seq_abc (int n_obs,
	      int n_part,
	      double *y1,
	      double *y2,
	      double *true_theta,
	      double **theta_part,
	      double tau,
	      bool noisy) {
  double toler = sqrt(tau);
  double *x1, *x2;
  double w[n_part];
  double prob_theta[n_part + 1];
  double *sim_y1[n_part];
  double *sim_y2[n_part];
  for (int i = 0; i < n_part; i++){
    sim_y1[i] = malloc(n_obs*sizeof(int));
    sim_y2[i] = malloc(n_obs*sizeof(int));
  }
  for (int j = 1; i < n_obs; j++){

    // Step 1
    for (int i = 0; i < n_part; i++)
      single_simul_obs(sim_y1[i],sim_y2[i],j,tau,theta_part[i],prob);
    
    // Step 2
    prob_theta[0] = 0;
    if (noisy == 0){
      gsl_rng_bivariate_gaussian(r,0.1,0.1,0,x1,x2);
      for (int i = 0; i < n_part; i++){
	w[i] = kernel(y1[j], y2[j], sim_y1[i][j] + toler * x1,
		      sim_y2[i][j] + toler * x2, toler);
	prob_theta[i+1] = prob_theta[i] + w[i];
      }
    } else {
      for (int i = 0; i < n_part; i++){
	w[i] = kernel(y1[j], y2[j], sim_y1[i][j], sim_y2[i][j], toler);
	prob_theta[i+1] = prob_theta[i] + w[i];
      }
    }
    
    // Step 3
    for (int i = 0; i < n_part; i++){
      int flag = -1;
      int l = 0;
      u = gsl_rng_uniform(r) * prob_theta[n_part];
      while (flag == -1){
	if (u <= prob[l+1] && u > prob[l]) flag = l;
	l++;
      }
      for (int k = 0; k < 3; k++)
	theta_part[i][k] = theta_part[l][k];
    }
  }
}

int main () {

  double true_theta[3] = {0.5, 0.0025, 0.3};
  int n_obs = 100;
  double tau = 0.1;
  double prob[3];
  int n_part = 5000;
  double h = sqrt(tau);
  int *y1 = malloc(sizeof(int)*n_obs);
  int *y2 = malloc(sizeof(int)*n_obs);

  double *theta_part[n_part];
  for (int i = 0; i < n_part; i++)
    theta_part[i] = malloc(sizeof(double)*3);

  y1[0] = 2;
  y2[0] = 2;  

  // Initialisazion
  simulation_obs(y1, y2, true_theta, prob, n_obs, h);
  prior_theta(theta_part, n_part);

  // Sequential ABC
  seq_abc(n_obs,y1,y2,true_theta,theta_part,tau,1);

  // Sequential noisy-ABC
  seq_abc(n_obs,y1,y2,true_theta,theta_part,tau,0);

  return 0;
}
