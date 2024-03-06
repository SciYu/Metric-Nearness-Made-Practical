#include "/usr/local/MATLAB/R2021a/extern/include/mex.h"
#include <string.h>
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    int n, niters, iter;
    int i, j, k;
    double *X, *X0, *D;
    double *Xi, *Xj, *Di;
    double v;

    /** get input **/
    X0 = (double*) mxGetData(prhs[0]);
    D = (double*) mxGetData(prhs[1]);
    niters = mxGetScalar(prhs[2]);
    n = mxGetM(prhs[0]);

    /** prepare output **/
    plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL); 
    X = (double*) mxGetData(plhs[0]);
    memcpy(X, X0, n*n*sizeof(double));

    for (iter = 0; iter < niters; iter++) {
        Xi = X; Di = D;
        for (i = 0; i < n; i++) {
            Xj = Xi;
            for (j=i+1; j<n; j++) {
                Xj += n;
                if (Xi[j] > Di[j]) {
                    Xi[j] = Di[j]; Xj[i] = Di[j];
                    for (k=0; k<n; k++) {
                        v = fabs(Xi[k]-Xj[k]);
                        if (v > Xi[j]) {
                            Xi[j] = v;
                        }
                    }
                } else {
                    Xi[j] = Di[j]; Xj[i] = Di[j];
                    for (k=0; k<n; k++) {
                        v = Xi[k] + Xj[k];
                        if (v < Xi[j]) {
                            Xi[j] = v;
                        }
                    }
                }
                Xj[i] = Xi[j];
            } // end of for j
            Xi += n; Di += n;
        } // end of for i
    } // end of for iter
}

