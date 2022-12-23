#include "functions.h"
#include<omp.h>
int maxi(int a,int b)
{
  if(a>b) return a;
  return b;
}

void check_matrix(double *mat, int n,int nrhs, int *flag)
{
    int i, j;
    double sum; 
    // flag == 2 for infinite solution
    // flag == 3 for No solution
    *flag = 2;
   // get to the last row 
   for (j=n;j<nrhs+n;j++){
    if(mat[(n-1)*(n+nrhs)+j]!=0) {*flag = 3; break;}
   }

}

/* function for exchanging two rows of
   a matrix */
void swap(double *mat, int row1, int row2,int n, int m)
{
  int max_dim = maxi(n,m);
  #pragma omp parallel for
    for (int i = 0; i < m; i++)
    {
        int temp = mat[row1*max_dim+i];
        mat[row1*max_dim+i] = mat[row2*max_dim+i];
        mat[row2*max_dim+i] = temp;
    }
}
void augmenter_matrix(double *mat ,double *a,double *b,int n,int nrhs)
{
    int i,j;
    int ldb = maxi(n,nrhs);
    #pragma omp parallel for private(j)
    for(i=0;i<n;i++){
      for(j=0;j<n+nrhs;j++){
        if(j<n) mat[i*(n+nrhs)+j]= a[i*n+j]; 
        else mat[i*(n+nrhs)+j] = b[i*ldb+(j-n)];
        }  
        }
   
}

//m = n+ nrhs
void Gauss_jordan(double *mat,int n ,int m, int *flag)
{
    int i,j,k,c;
    int max_dim = maxi(n,m); 
    double norm;
    for(i=0;i<n;i++){
        if(mat[i*(max_dim+1)]==0){ // if mat[i][i] = 0 search a row l (l>i) where mat[l][l]!=0
          c = 1;
            while ((i + c) < n && mat[(i+c)*max_dim+i] == 0)
                c++;        
            if ((i + c) == n) {
                *flag = 1;
                break;
            }
            swap(mat,i,i+c,n,m); // swap the two rows
            }
        for(j=0;j<n;j++){
          if(i!=j){
            norm = mat[j*max_dim+i]/mat[i*(max_dim+1)];
            //#pragma omp parallel for
            for(k=0;k<m;k++){
                mat[j*max_dim+k]-=norm*mat[i*max_dim+k];
            }
      
        }
    }
}
}
