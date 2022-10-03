#include <time.h>
#include<unistd.h>
# include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>
#include "functions.h"
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


int my_dgesv(int n, int nrhs, double *a, double *b)
{
      int ldb = maxi(n,nrhs);
      double *mat = (double *)calloc(sizeof(double),n*(n+nrhs));
      augmenter_matrix(mat,a,b,n,nrhs);
      int i,j,flag=0;
      Gauss_jordan(mat,n,n+nrhs,&flag);
      if(flag==0){ // there is an unique solution
         for(i=0;i<n;i++){
          for(j=0;j<nrhs;j++){
              b[i*ldb+j]=mat[i*(n+nrhs)+(n+j)]/mat[i*((n+nrhs)+1)];
          }
        }
      } else{
        check_matrix(mat,n,nrhs,&flag);   
      }   
free(mat);
return flag;
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

  info = my_dgesv(n, nrhs, a, b);
  printf("Time taken by my implementation: %.2fs\n", (double) (clock() - tStart) / CLOCKS_PER_SEC);
  if(info==0){
      if (check_result(bref, b, size) == 1)
        printf("Result is ok!\n");
      else
        printf("Result is wrong!\n");
      }else if (info==2)
      {
        printf("Infinite solutions\n");
      }else printf("No solution\n");
      
  free(ipiv);free(ipiv2);
}
