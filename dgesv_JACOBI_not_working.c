#include <time.h>
#include<unistd.h>
# include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>
// <mkl_lapacke.h>

double *generate_matrix(int size)
{
  int i;
  double *matrix = (double *) malloc(sizeof(double) * size * size);

  srand(1);

  for (i = 0; i < size * size; i++) {
    matrix[i] = rand() % 100;
  }

  return matrix;
}

int is_nearly_equal(double x, double y)
{
  const double epsilon = 1e-5 /* some small number */;
  return abs(x - y) <= epsilon * abs(x);
  // see Knuth section 4.2.2 pages 217-218
}

int check_result(double *bref, double *b, int size)
{
  int i;

  for(i = 0; i < size*size; i++) {
    if (!is_nearly_equal(bref[i], b[i]))
      return 0;
  }

  return 1;
}

int check(double *bref, double *b, int size)
{
  int i;

  for(i = 0; i < size; i++) {
    if (fabs(bref[i]-b[i])>pow(10,-2))
      return 0;
  }

  return 1;
}
double norme(double *u,int size){
    double sum = 0.0;
    int i;
    for(i=0;i<size;i++)
      sum+=u[i]*u[i];
    return sqrt(sum);
}
int my_dgesv(int n, int nrhs, double *A, int lda, int *ipiv, double *b, int ldb)
{

  //Replace next line to use your own DGESV implementation
  //LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, a, lda, ipiv, b, ldb);
  int boolean ;
  double *D = (double *)calloc(sizeof(double),n*n);
  double sum;
  double *N = (double *)calloc(sizeof(double),n*n);
  double *x = (double *)calloc(sizeof(double),n*nrhs);
  double *tmp = (double *)calloc(sizeof(double),n*nrhs); // N*x +b n*nhrs
  double *tmp2 = (double *)calloc(sizeof(double),n*nrhs); // y-x n*nhrs
  double *y = (double *)calloc(sizeof(double),n*nrhs); // D-1(N*x +b)
  int i,j,k;
  // compute D
  for(i=0;i<n;i++) D[i*n+i]=A[i*n+i]; 
// compute N
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      if(i!=j)
        N[i*n+j]=-A[i*n+j];
     }
  }



// initial X0
  for(i=0;i<n;i++)
    for(j=0;j<nrhs;j++)
       x[i*n+j] = (double)rand()/(double) RAND_MAX;


  do{
    boolean = 0;
    for(i=0;i<n;i++){
      for(j=0;j<nrhs;j++){
        tmp[i*n+j]=0;
        for(k=0;k<n;k++)
          tmp[i*n+j] += N[i*n+k]*x[k*n+j];
        }
    }
    for(i=0;i<n;i++){
          for(j=0;j<nrhs;j++){
            tmp[i*n+j]+=b[i*n+j];
            }
        }

//compute D-1*(Nx+b)
    for(i=0;i<n;i++){
      for(j=0;j<nrhs;j++){
        y[i*n+j] = 0;
        for(k=0;k<n;k++){
          if(D[i*n+k]!=0)
          y[i*n+j] += tmp[k*n+j]/D[i*n+k];}
        }
      }
    printf("y %f\t",y[0]);
  printf("y %f\t",y[1]);
  printf("y %f\t",y[2]);
  printf("y %f\n",y[3]);
    printf("x %f\t",x[0]);
  printf("x %f\t",x[1]);
  printf("x %f\t",x[2]);
  printf("x %f\n",x[3]);


  boolean = check(y,x,n*nrhs);
  sleep(2);
    for(i=0;i<n;i++)
    for(j=0;j<nrhs;j++)
      x[i*n+j] = y[i*n+j];
  }while(!boolean);

  for(i=0;i<n;i++)
    for(j=0;j<nrhs;j++)
      b[i*n+j] = y[i*n+j];

free(D);free(N);free(tmp);free(tmp2);
free(x);free(y);
}

void main(int argc, char *argv[])
{
  int size = atoi(argv[1]);

  double *a, *aref;
  double *b, *bref;

  a = generate_matrix(size);
  aref = generate_matrix(size);
  b = generate_matrix(size);
  bref = generate_matrix(size);
  // Using LAPACK dgesv OpenBLAS implementation to solve the system
  int n = size, nrhs = size, lda = size, ldb = size, info;
  int *ipiv = (int *) malloc(sizeof(int) * size);

  clock_t tStart = clock();
  info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, aref, lda, ipiv, bref, ldb);
  printf("Time taken by OpenBLAS LAPACK: %.2fs\n", (double) (clock() - tStart) / CLOCKS_PER_SEC);

  int *ipiv2 = (int *) malloc(sizeof(int) * size);

  tStart = clock();
  my_dgesv(n, nrhs, a, lda, ipiv2, b, ldb);
  /*printf("%f\t",b[0]);
  printf("%f\t",b[1]);
  printf("%f\t",b[2]);
  printf("%f\t",b[3]);*/
  printf("Time taken by my implementation: %.2fs\n", (double) (clock() - tStart) / CLOCKS_PER_SEC);
  if (check_result(bref, b, size) == 1)
    printf("Result is ok!\n");
  else
    printf("Result is wrong!\n");

}
