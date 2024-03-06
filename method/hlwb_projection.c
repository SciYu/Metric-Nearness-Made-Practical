#include "/usr/local/MATLAB/R2021a/extern/include/mex.h"
#include <string.h>
#include "time.h"

// X = lambda * D + (l-lambda) * X 
void addm(double *D, double *X, double lambda, int n) {
    double lambda2 = 1.0 - lambda;
    for (int i=0; i<n*n; i++) {
        X[i] = lambda * D[i] + lambda2 * X[i];
    }
}

// nmse between two matrices, i.e., ||X-D||^2 / ||X||^2
double compute_nmse(double *D, double *X, int n) {
    double se = 0.0;
    double normD = 0.0;
    for (int i=0; i<n*n; i++) {
        double tmp = D[i] - X[i];
        se += (tmp * tmp);
        normD += D[i]*D[i]; 
    }
    return se / normD;
}

// number of vltns to triangle inequality
double compute_vltns(double *X, int n) {
    double ab;
    double count = 0.0;
    double* Xi = X;
    const double epsilon = 1.0e-10;
    for (int i=0; i<n; i++) {
        double* Xj = Xi;
        for (int j=i+1; j<n; j++) {
            Xj += n;
            ab = Xi[j];
            for (int k=0; k<n; k++) {
                if (k==i || k==j) continue;
                if (ab - epsilon > Xj[k] + Xi[k]) {
                    count += 1.0;
                }
            }
        }
        Xi += n;
    }
    return count;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    int n, niters, iter;
    int i, j, k;
    double *X, *X0, *D0;
    double *Xi, *Xj, *Xk;
    double delta;
    const double epsilon = 1e-10;

    /** get input **/
    X0 = (double*) mxGetData(prhs[0]);
    D0 = (double*) mxGetData(prhs[1]);
    niters = mxGetScalar(prhs[2]);
    n = mxGetM(prhs[0]);

    /** prepare output **/
    plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL); 
    X = (double*) mxGetData(plhs[0]);

    /** print info **/
    memcpy(X, X0, n*n*sizeof(double));
    clock_t elapse = clock();
    double updates = 0.0d;
    mexPrintf("Data matrix. %d x %d\n", n, n);
    mexPrintf("HLWB projection Entered\n");
    mexPrintf("Iters\tNMSE_before\tVLTN_before\tUpdates  \tNMSE_after\tVLTN_after\tSecs\t\n");

    for (iter=1; iter<=niters; iter++) {
        addm(D0, X, 0.382/iter, n);
        mexPrintf("%d\t%0.7f\t%0.3fK  \t%0.6fM\t", iter, compute_nmse(D0,X,n), compute_vltns(X,n)/1.0e+3, updates/1.0e+6);
        clock_t elapse = clock();
        Xi = X;
        for (i=0; i<n; i++) {
            Xj = Xi;
            for (j=i+1; j<n; j++) {
                Xj += n; Xk = Xj;
                for (k=j+1; k<n; k++) {
                    Xk += n;
                    delta = Xi[j]-Xi[k]-Xj[k];
                    if (delta > 0.0) {
                        delta /= 3.0;
                        Xi[j] -= delta; Xj[i]=Xi[j];
                        Xi[k] += delta; Xk[i]=Xi[k];
                        Xj[k] += delta; Xk[j]=Xj[k];
                        updates += 3.0;
                        continue;
                    }
                    delta = Xi[k]-Xi[j]-Xj[k];
                    if (delta > 0.0) {
                        delta /= 3.0;
                        Xi[k] -= delta; Xk[i]=Xi[k];
                        Xi[j] += delta; Xj[i]=Xi[j];
                        Xj[k] += delta; Xk[j]=Xj[k];
                        updates += 3.0;
                        continue;
                    }
                    delta = Xj[k]-Xi[j]-Xi[k];
                    if (delta > 0.0) {
                        delta /= 3.0;
                        Xj[k] -= delta; Xk[j]=Xj[k];
                        Xi[j] += delta; Xj[i]=Xi[j];
                        Xi[k] += delta; Xk[i]=Xi[k];
                        updates += 3.0;
                        continue;
                    }
                }
            } // end of for j
            Xi += n;
        } // end of for i
        double secs = (clock()-elapse)*1.0/CLOCKS_PER_SEC;
        mexPrintf("%0.7f\t%0.3fK  \t%0.3f\n", compute_nmse(D0,X,n), compute_vltns(X,n)/1.0e+3, secs);
    }

    // zero diagonal and make symmetric
    Xi = X;
    for (i=0; i<n; i++) {
        Xi[i] = 0.0d;
        Xj = Xi;
        for (j=i+1; j<n; j++) {
            Xj += n;
            Xj[i] = Xi[j];
        }
        Xi += n;
    }
}

