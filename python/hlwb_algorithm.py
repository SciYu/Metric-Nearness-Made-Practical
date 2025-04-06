import numpy as np
import OptimizationPack.matrix_optimization as matrix_optimization
# from numba import njit, prange
# from scipy.linalg import eigh

def nearpsd(A, maxits=100, low=0, high=1, d=1):
    """
    Computes the nearest positive semi-definite matrix for a given square matrix.

    Parameters:
    A : numpy.ndarray
        A square matrix to be calibrated.
    maxits : int, optional
        Maximum number of iterations allowed, default is 100.
    low : float, optional
        Minimum value for elements, default is 0.
    high : float, optional
        Maximum value for elements, default is 1.
    d : float, optional
        Value of the diagonal elements, default is 1.

    Returns:
    numpy.ndarray
        Nearest positive semi-definite matrix to A.
    """
    if not np.allclose(A, A.T):
        A = (A + A.T) / 2
    n = A.shape[0]
    
    U = np.zeros((n, n))
    Y = A.copy()
    tolconv = 1.0e-6
    toleigs = 1.0e-5
    iter = 0

    while True:
        T = Y - U

        # Project onto PSD matrices
        eigvals, eigvecs = np.linalg.eigh(T)
        p = eigvals > toleigs
        X = eigvecs[:, p] @ np.diag(eigvals[p]) @ eigvecs[:, p].T

        # Update correction
        U = X - T

        # Maximum iteration and convergence test
        iter += 1
        if iter >= maxits or np.linalg.norm(Y - X, np.inf) / np.linalg.norm(Y, np.inf) <= tolconv:
            break

        # Update Y with constraints
        Y = X
        np.fill_diagonal(Y, d)
        Y = np.clip(Y, low, high)

    np.fill_diagonal(Y, d)
    Y = np.clip(Y, low, high)

    return Y


def project_to_kernel(D, mu=1, maxits=100, tolconv=1.0e-6, toleigs=1.0e-5):
    """
    Projects a distance matrix to a new distance matrix with a positive semi-definite (PSD) kernel matrix.
    
    Parameters:
    D : numpy.ndarray
        Pairwise distance matrix to be calibrated.
    mu : float, optional
        Calibration parameter, default is 1.
    maxits : int, optional
        Maximum number of iterations for PSD projection, default is 100.
    tolconv : float, optional
        Convergence tolerance, default is 1.0e-5.
    toleigs : float, optional
        Eigenvalue tolerance for PSD projection, default is 1.0e-5.
    
    Returns:
    numpy.ndarray
        Calibrated distance matrix that induces a positive semi-definite kernel.
    """
    
    # Transform distance matrix to the kernel matrix
    gamma = -mu / np.max(D)
    K = np.exp(gamma * D)
    n = K.shape[0]
    low_val = np.exp(-mu)
    high_val = 1
    diag_val = 1

    Y = (K + K.T) / 2
    U = np.zeros((n, n))
    iter_count = 0

    # Iteratively project onto the nearest PSD matrix
    while True:
        T = Y - U
        
        # Project onto PSD matrices by keeping only positive eigenvalues
        eigvals, eigvecs = np.linalg.eigh(T)
        pos_eigvals = eigvals > toleigs
        X = eigvecs[:, pos_eigvals] @ np.diag(eigvals[pos_eigvals]) @ eigvecs[:, pos_eigvals].T

        # Update correction and check for convergence
        U = X - T
        iter_count += 1
        if iter_count >= maxits or np.linalg.norm(Y - X, np.inf) / np.linalg.norm(Y, np.inf) <= tolconv:
            break

        # Update Y with constraints, maintaining symmetry and clipping
        Y = X
        np.fill_diagonal(Y, diag_val)
        Y = np.clip(Y, low_val, high_val)

    # Step 5: Convert back to the distance matrix
    C = np.log(Y) / gamma
    np.fill_diagonal(C, 0) # Ensure zero-diagonal
    D_new = (C + C.T) / 2  # Ensure symmetry

    print(iter_count)
    
    return D_new


def hlwb_algorithm(D, n_projection=100):
    """
    The alternating projection stage iteratively refines the approximate solution X0 to the optimal X.

    Parameters:
    D : numpy.ndarray
        Initial non-metric distance matrix.
    n_projection : int, optional
        Number of projections, default is 100.

    Returns:
    numpy.ndarray
        Optimal solution to the metric nearness problem.
    """

    D_proj = project_to_kernel(D, maxits=10)
    
    X0 = matrix_optimization.heuristic_improve(D_proj.copy(), D.copy(), n_improve=1)

    X = matrix_optimization.hlwb_projection(X0.copy(), D.copy(), n_projection=n_projection)
    
    X = (X + X.T) / 2
    np.fill_diagonal(X, 0)

    return X