#ifndef __functions_H_
#define __functions_H_
#include "functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/**
* \fn maxi(int a,int b)
* \author Manuel NKEGOUM<nkegoumnzo@cy-tech.fr>
*\version 0.1
*\date 2/10/2022
*\param a integer
*\param b integer
 find the maximum between two integers
*/
int maxi(int a,int b);


/**
* \fn check_matrix(double *mat, int n,int nrhs, int *flag);
* \author Manuel NKEGOUM<nkegoumnzo@cy-tech.fr>
*\version 0.1
*\date 2/10/2022
*\param mat augmented matrix
*\param n number of rows of the augmented matrix mat
*\param nrhs augmented matrix has n+nrhs columns
*\param flag represent the consistency (2 => infinite solutions, 3 => No solution)
 check the consistency of the system of equations
*/
void check_matrix(double *mat, int n,int nrhs, int *flag);


/**
* \fn swap(double *mat, int row1, int row2,int n, int m)
* \author Manuel NKEGOUM<nkegoumnzo@cy-tech.fr>
*\version 0.1
*\date 2/10/2022
*\param mat augmented matrix
*\param row1 row1 to swap 
*\param row2 row2 to swap 
*\param n number of rows of the augmented matrix mat
*\param m number of columns of the augmented matrix mat
 swap two rows in a matrix
*/
void swap(double *mat, int row1, int row2,int n, int m);

/**
* \fn augmenter_matrix(double *mat ,double *a,double *b,int n,int nrhs)
* \author Manuel NKEGOUM<nkegoumnzo@cy-tech.fr>
*\version 0.1
*\date 2/10/2022
*\param mat augmented matrix [a|b]
*\param a matrix A 
*\param b matrix B
*\param n number of rows of A
*\param nrhs number of columns of B
initialize the augmented matrix which should be n x (n+nrhs)
*/
void augmenter_matrix(double *mat ,double *a,double *b,int n,int nrhs);


/**
* \fn Gauss_jordan(double *mat,int n ,int m, int *flag)
* \author Manuel NKEGOUM<nkegoumnzo@eisti.eu>
*\version 0.1
*\date 10 décembre 2020
*\param mat augmented matrix
*\param n number of rows of the augmented matrix mat
*\param m number of columns of the augmented matrix mat
*\param flag represent the consistency (2 => infinite solutions, 3 => No solution)
Apply the Reduction of Gauss-Jordan upon the augmented matrix
*/
void Gauss_jordan(double *mat,int n ,int m, int *flag);

#endif