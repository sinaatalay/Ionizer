#pragma once

// FROM: https://people.sc.fsu.edu/~jburkardt/cpp_src/mgmres/mgmres.html

void atx_st(int n, int nz_num, int ia[], int ja[], double a[], double x[], double w[]);
void ax_st(int n, int nz_num, int ia[], int ja[], double a[], double x[], double w[]);

void mgmres_st(int n, int nz_num, int ia[], int ja[], double a[], double x[], double rhs[], int itr_max, int mr, double tol_abs, double tol_rel);
void mult_givens(double c, double s, int k, double g[]);

double r8vec_dot(int n, double a1[], double a2[]);
double* r8vec_uniform_01(int n, int& seed);
