import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg

class WgMode:
    def __init__(self, beta, field):
        self.beta = beta
        self.field = field

def solve(omega, epsilon, order=1):
    """ Form the matrices we will need to solve the eigenvalue problem.

    """
    n = epsilon.size 
    diag_vals = [[1], [-1], [-1]] * np.ones((3, n))
    Dx = sparse.spdiags(diag_vals, (0, -1, n-1), n, n)

    epsilon_x = epsilon
    epsilon_y = 0.5 * (epsilon + epsilon[[n-1] + range(n-1)])

    my_diag = lambda x: sparse.spdiags(x, 0, n, n)

    A1 = -1 * my_diag(epsilon_x) * Dx.T * my_diag(epsilon_y**-1) * Dx;
    A2 = my_diag(epsilon_x)

    A1 = A1 + 4 * sparse.eye(n, n)

    # Calculate relevant eigenvalues.
    beta2, fields= linalg.eigs(A1 + omega**2 * A2, k=order)
    
    # Find the mode with the lowest beta2, which is the one we want.
    ind = np.argmin(beta2)

    return WgMode(np.sqrt(beta2[ind] - 4.0), fields[:, ind])
