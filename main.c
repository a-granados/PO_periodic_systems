#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <math.h>
#include "rk78.h" //Runge-Kutta integrator
#include "memory.h" //Memory allocation for matrices and vectors
#include "lu.h" //LU method
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>


using std::cout;
using std::endl;
using namespace std;

//Global variables:
double twopi=8*atan(1.);
int n=2; //Dimension of the system
double T=twopi;//Period of the non-autonomous system
double t0=0;



void vfield(double t, double *x, int ndim, double *dx);
int Newton_method (int q,double *xout,double **DPout);
int strobomap(double* xini, double tf, double **DP);
void stability(double **A);

int main(int argc, char * argv[]){
  int i;
  double *x=new double[n];
  double **DP=new double*[n];
  for (i=0;i<n;i++){
    DP[i]=new double[n];
  }
  if (argc!=n+1){
    cout <<"Wrong use, provide "<<n<<" initial conditions"<<endl;
  }
  for (i=0;i<n;i++){
    sscanf(argv[i+1], "%lg", &x[i]);
  }

  Newton_method(1,x,DP);

  
  cout <<"Periodic point:"<<endl;
  for (i=0;i<n;i++){
    cout<<x[i]<<" ";
  }
  cout <<endl;
  cout <<"Eigenvalues and eigenvectors of Ds"<<endl;
  stability(DP);

  delete[] x;
  for (i=0;i<n;i++){
    delete[] DP[i];
  }
  delete[] DP;
}


int Newton_method(int q,double *xout,double **DPout){

  //This routine peforms a Newton method to compute a q-periodic point of
  //the stroboscopic map. It returns 0 if the method converged, and 1
  //otherwise.
  //
  //Input parameters:
  //- p: period of the periodic point of the stroboscopic map that we
  //are looking for. For example, if q=1 we are then looking for a
  //fixed point of the stroboscopic map; that is, an initial condition
  //for a T-periodic orbit of the original system.
  //- xout: contains an initial guess where to start the Newton method.
  //
  //Output parameters:
  //- xout. If  the method converged, then the n-periodic point is stored in xout.
  //- DPout: If the method converged, the routine returns the
  //differential of the stroboscopic map evaluated at the fixed point.
  //This is useful to compute stability properties, eigenvalues and
  //eigenvectors.

  double tf;
  int i,j;
  double *xa;
  int *perm;
  double tol=1.0e-12;
  double NewtonTol=1.0e-8;
  double dist;
  int nite=0;
  int result;
  double *p;
  double *dta;
  int maxiterNewton=500;
  int res=0;
  double *x=new double[n];
  double **DP;
  DP=crear_matriu(0,n-1,0,n-1);
  for (i=0;i<n;i++){
    x[i]=xout[i];
  }

  p=(double *)calloc(n, sizeof(double));
  dta=(double *)calloc(n, sizeof(double));
  xa=(double *)calloc(n, sizeof(double));
  perm=(int *)calloc(n, sizeof(int));;

  tf=double(q)*T;

  for (i=0;i<n;i++){
  xa[i]=x[i];
  }  
/*integrate*/
  
  res=strobomap(x, tf, DP);
  if (res==1){
    free(xa);
    free(perm);
    free(dta);
    free(p);
    alliberar_matriu(DP,0,n-1,0);
    delete[] x;
    return 1;
  }
  dist=0;
  for (i=0;i<n;i++){
    dist+=(x[i]-xa[i])*(x[i]-xa[i]);
  }
  dist=sqrt(dist);
  
  while(dist > NewtonTol && nite<maxiterNewton){

     for(i=0; i<n; i++){
        perm[i]=i;
     }
   
     for(i=0;i<n;i++){
        DP[i][i]=DP[i][i]-1.;
        p[i]=-(x[i]-xa[i]);
     }

     result=lu(DP,n,perm,tol);
     if (result==0){
       cout <<"# LU failed."<<endl;
	//exit(1);
	free(xa);
	free(perm);
	free(dta);
	free(p);
	alliberar_matriu(DP,0,n-1,0);
	delete[] x;
	return 1;
     }
     
     solve(DP, dta, p, n, perm);
     
     for (i=0;i<n;i++){
       xa[i]=xa[i]+dta[i];
     }
     
     for (i=0; i<n; i++){
        x[i]=xa[i];
     }

     res=strobomap(x, tf, DP);
    if (res==1){
      free(xa);
      free(perm);
      free(dta);
      free(p);
      alliberar_matriu(DP,0,n-1,0);
      delete[] x;
      return 1;
    }

     dist=0;
     for (i=0;i<n;i++){
       dist+=(x[i]-xa[i])*(x[i]-xa[i]);
     }
     dist=sqrt(dist);

     nite++;
  }
  free(xa);
  free(perm);
  free(dta);
  free(p);

   if (nite==maxiterNewton && dist > NewtonTol || isnan(x[0]) ||isnan(x[1])){
     //printf("#massa iteracions newton\n");
     //exit(1);
     alliberar_matriu(DP,0,n-1,0);
     delete[] x;
     return 1;
   }
   else {
     for (i=0;i<n;i++){
       xout[i]=x[i];
       for (j=0;j<n;j++){
	 DPout[i][j]=DP[i][j];
       }
     }
     alliberar_matriu(DP,0,n-1,0);
     delete[] x;
     return 0;
   }
}

int strobomap(double* xini, double tf, double **DP){
  double t=t0; 
  int i,j;
  int ndim=n+n*n;//Dimension of the system with the variational equations
  double *dx=new double[ndim];
  double *x=new double[ndim];
  double hini=0.05;
  double h=hini;
  double hmin=1.0E-4;
  double hmax=1.0;
  double tol=1.0e-13;

  for (i=0;i<n;i++){
    x[i]=xini[i];
  }
  //Identity matrix is the initial condition for the variational
  //equations.
  for (i=n;i<ndim;i++){
    if ((i-n)%(n+1)==0){
      x[i]=1.;
    }
    else{
      x[i]=0.;
    }
  }
 ini_rk78(ndim);
 while (t<t0+tf){
  rk78(&t,x,&h,tol,hmin,hmax,ndim,vfield);
 }

 h=-(t-tf-t0);
 rk78(&t,x,&h,tol,fabs(h),fabs(h),ndim,vfield);
 end_rk78(ndim);

 for (i=0;i<n;i++){
   xini[i]=x[i];
 }

 for (i=0;i<n;i++){
   for (j=0;j<n;j++){
     DP[i][j]=x[(i+1)*n+j];
   }
 }

 delete[] dx;
 delete[] x;
 return 0;
}

void vfield(double t, double *x, int ndim, double *dx){
  //Equations of the pendulum:
  double eps=5e-2;
  dx[0]=x[1];
  dx[1]=-sin(x[0])+eps*sin(twopi/T*t);

  //Variational equations:
  dx[2]=x[4];
  dx[3]=x[5];
  dx[4]=-cos(x[0])*x[2];
  dx[5]=-cos(x[0])*x[3];
}

void stability(double **A){
  //This computes and prints the eigenvalues and eigenvectors of the
  //matrix A
  gsl_matrix *M=gsl_matrix_alloc(n,n);
  gsl_vector_complex *eval = gsl_vector_complex_alloc (n);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc (n, n);
  int i,j;

  for (i=0;i<n;i++){
    for (j=0;j<n;j++){
      gsl_matrix_set(M,i,j,A[i][j]);
    }
  }

  gsl_eigen_nonsymmv_workspace * w = 
    gsl_eigen_nonsymmv_alloc (n);
  
  gsl_eigen_nonsymmv (M, eval, evec, w);

  gsl_eigen_nonsymmv_free (w);
  gsl_matrix_free(M);

  gsl_eigen_nonsymmv_sort (eval, evec, 
                           GSL_EIGEN_SORT_ABS_ASC);
  
  {

    for (i = 0; i < n; i++)
      {
        gsl_complex eval_i 
           = gsl_vector_complex_get (eval, i);
        gsl_vector_complex_view evec_i 
           = gsl_matrix_complex_column (evec, i);
        printf ("eigenvalue = %g + %gi\n",
         GSL_REAL(eval_i), GSL_IMAG(eval_i));
        printf ("eigenvector = \n");
	//cout <<GSL_REAL(eval_i)<<" + "<<GSL_IMAG(eval_i)<<" i"<<endl;
        for (j = 0; j < n; ++j)
          {
            gsl_complex z = 
              gsl_vector_complex_get(&evec_i.vector, j);
         printf("%g + %gi\n", GSL_REAL(z), GSL_IMAG(z));
          }
      }
  }
}

