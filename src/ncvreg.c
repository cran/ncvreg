#include <math.h>
#include <string.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <R_ext/Applic.h>

static double *vector(int n)
{
  double *v;
  v = Calloc(n, double);
  return v;
}

static void free_vector(double *v)
{
  Free(v);
}

static int checkConvergence(double *beta, double *beta_old, double eps, int l, int J)
{
  int j;
  int converged = 1;
  for (j=0; j < J; j++)
    {
      if (beta[l*J+j]!=0 & beta_old[j]!=0)
	{
	  if (fabs((beta[l*J+j]-beta_old[j])/beta_old[j]) > eps)
	    {
	      converged = 0;
	      break;
	    }
	}
      else if (beta[l*J+j]==0 & beta_old[j]!=0)
	{
	  converged = 0;
	  break;
	}
      else if (beta[l*J+j]!=0 & beta_old[j]==0)
	{
	  converged = 0;
	  break;
	}
    }
  return(converged);
}

static double S(double z, double l)
{
  if (z > l) return(z-l);
  if (z < -l) return(z+l);
  return(0);
}

static double MCP(double z, double l, double a, double v)
{
  double s=0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  if (fabs(z) <= l) return(0);
  if (fabs(z) <= a*l) return(s*(fabs(z)-l)/(v*(1-1/a)));
  return(z/v);
}

static double SCAD(double z, double l, double a, double v)
{
  double s=0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  if (fabs(z) <= l) return(0);
  if (fabs(z) <= 2*l) return(s*(fabs(z)-l)/v);
  if (fabs(z) <= a*l) return(s*(fabs(z)-a*l/(a-1))/(v*(1-1/(a-1))));
  return(z/v);
}

static void cdfit_gaussian(double *beta, int *iter, double *x, double *y, int *n_, int *p_, char **penalty_, double *lambda, int *L_, double *eps_, int *max_iter_, double *a_)
{
  /* Declarations */
  int L=L_[0]; int p=p_[0]; int n=n_[0]; int max_iter=max_iter_[0]; double eps=eps_[0]; double a=a_[0]; char *penalty=penalty_[0];
  int converged;
  double *r, *beta_old;
  r = vector(n); for (int i=0;i<n;i++) r[i] = y[i];
  beta_old = vector(p);

  /* Path */
  for (int l=1;l<L;l++)
    {
      for (int j=0;j<p;j++) beta_old[j] = beta[(l-1)*p+j];
      while (iter[l] < max_iter)
	{
	  converged = 0;
	  iter[l] = iter[l] + 1;

	  /* Covariates */
	  for (int j=0;j<p;j++)
	    {
	      /* Calculate z */
	      double z = 0;
	      for (int i=0;i<n;i++) z = z + x[j*n+i]*r[i];
	      z = z/n + beta_old[j];

	      /* Update beta_j */
	      if (strcmp(penalty,"MCP")==0) beta[l*p+j] = MCP(z,lambda[l],a,1);
	      if (strcmp(penalty,"SCAD")==0) beta[l*p+j] = SCAD(z,lambda[l],a,1);

	      /* Update r */
	      if (beta[l*p+j] != beta_old[j]) for (int i=0;i<n;i++) r[i] = r[i] - (beta[l*p+j] - beta_old[j])*x[j*n+i];
	    }

	  /* Check for convergence */
	  if (checkConvergence(beta,beta_old,eps,l,p))
	    {
	      converged  = 1;
	      break;
	    }
	  for (int j=0;j<p;j++) beta_old[j] = beta[l*p+j];
	}
      /*if (converged==0) warning("Failed to converge");*/
    }

  free_vector(beta_old);
  free_vector(r);
}

static void cdfit_binomial(double *beta0, double *beta, int *iter, double *x, double *y, int *n_, int *p_, char **penalty_, double *lambda, int *L_, double *eps_, int *max_iter_, double *a_)
{
  /* Declarations */
  int L=L_[0];int p=p_[0];int n=n_[0];int max_iter=max_iter_[0];double eps=eps_[0];double a=a_[0];char *penalty=penalty_[0];
  int converged;
  double beta0_old;
  double *r, *w, *beta_old;
  r = vector(n);
  w = vector(n);
  beta_old = vector(p);

  /* Initialization */
  double ybar=0;
  for (int i=0;i<n;i++) ybar = ybar + y[i];
  ybar = ybar/n;
  beta0[0] = log(ybar/(1-ybar));

  /* Path */
  double xwr,xwx,eta,pi,yp,yy,z,v;
  for (int l=1;l<L;l++)
    {
      beta0_old = beta0[l-1];
      for (int j=0;j<p;j++) beta_old[j] = beta[(l-1)*p+j];
      while (iter[l] < max_iter)
	{
	  converged = 0;
	  iter[l] = iter[l] + 1;

	  /* Approximate L */
	  yp=yy=0;
	  for (int i=0;i<n;i++)
	    {
	      eta = beta0_old;
	      for (int j=0;j<p;j++) eta = eta + x[j*n+i]*beta_old[j];
	      pi = exp(eta)/(1+exp(eta));
	      if (pi > .9999)
		{
		  pi = 1;
		  w[i] = .0001;
		}
	      else if (pi < .0001)
		{
		  pi = 0;
		  w[i] = .0001;
		}
	      else w[i] = pi*(1-pi);
	      r[i] = (y[i] - pi)/w[i];
	      yp = yp + pow(y[i]-pi,2);
	      yy = yy + pow(y[i]-ybar,2);
	    }
	  if (yp/yy < .01)
	    {
	      warning("Model saturated; exiting...");
	      for (int ll=l;ll<L;ll++)
		{
		  beta0[ll] = R_NaReal;
		  for (int j=0;j<p;j++) beta[ll*p+j] = R_NaReal;
		}
	      free_vector(beta_old);
	      free_vector(w);
	      free_vector(r);
	      return;
	    }

	  /* Intercept */
	  xwr = xwx = 0;
	  for (int i=0;i<n;i++)
	    {
	      xwr = xwr + w[i]*r[i];
	      xwx = xwx + w[i];
	    }
	  beta0[l] = xwr/xwx + beta0_old;
	  for (int i=0;i<n;i++) r[i] = r[i] - (beta0[l] - beta0_old);

	  /* Covariates */
	  for (int j=0;j<p;j++)
	    {
	      /* Calculate z */
	      xwr=0;
	      xwx=0;
	      for (int i=0;i<n;i++)
		{
		  xwr = xwr + x[j*n+i]*w[i]*r[i];
		  xwx = xwx + x[j*n+i]*w[i]*x[j*n+i];
		}
	      z = xwr/n + (xwx/n)*beta_old[j];
	      v = xwx/n;

	      /* Update beta_j */
	      if (strcmp(penalty,"MCP")==0) beta[l*p+j] = MCP(z,lambda[l],a,v);
	      if (strcmp(penalty,"SCAD")==0) beta[l*p+j] = SCAD(z,lambda[l],a,v);

	      /* Update r */
	      if (beta[l*p+j] != beta_old[j]) for (int i=0;i<n;i++) r[i] = r[i] - (beta[l*p+j] - beta_old[j])*x[j*n+i];
	    }

	  /* Check for convergence */
	  if (checkConvergence(beta,beta_old,eps,l,p))
	    {
	      converged = 1;
	      break;
	    }
	  beta0_old = beta0[l];
	  for (int j=0;j<p;j++) beta_old[j] = beta[l*p+j];
	}
    }

  free_vector(beta_old);
  free_vector(w);
  free_vector(r);
}

static const R_CMethodDef cMethods[] = {
  {"cdfit_gaussian", (DL_FUNC) &cdfit_gaussian, 12},
  {"cdfit_binomial", (DL_FUNC) &cdfit_binomial, 13},
  NULL
};

void R_init_ncvreg(DllInfo *info)
{
  R_registerRoutines(info,cMethods,NULL,NULL,NULL);
}
