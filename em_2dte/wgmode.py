import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg

class WgMode:
    """ Class for working with waveguide modes.

    """

    def __init__(self, beta, field):

        self.beta = beta
        self.field = field

        # If beta is real (i.e. we have a true propagating mode), 
        # then make sure that the fields are all real as well.
        if np.isreal(self.beta):
            if all(np.isreal(self.field)):
                self.field = np.real(self.field)
            else:
                raise Exception("Waveguide mode should be real-valued.")

    def plot(self, gplot):
        gplot.plot(self.field)

def solve(omega, epsilon, order=1):
    """ Form the matrices we will need to solve the eigenvalue problem.

    """
    n = epsilon.size # Grid size (1D).

    # Create the derivative matrix (periodic boundary conditions).
    diag_vals = [[1], [-1], [-1]] * np.ones((3, n))
    Dx = sparse.spdiags(diag_vals, (0, -1, n-1), n, n)

    # Calculate the permittivity, assuming epsilon is located at (0, 0, 0) on 
    # the Yee grid.
    epsilon_x = epsilon
    epsilon_y = 0.5 * (epsilon + epsilon[[n-1] + range(n-1)])

    # Helper function to quickly create a sparse, diagonal matrix.
    my_diag = lambda x: sparse.spdiags(x, 0, n, n)

    # Form A. Shift the matrix by 4 so that all eigenvalues are positive.
    A = -1 * my_diag(epsilon_x) * Dx.T * my_diag(epsilon_y**-1) * Dx + \
        omega**2 * my_diag(epsilon_x) + \
        4 * sparse.eye(n, n)

    # Calculate relevant eigenvalues.
    # These will be the maximal eigenvalues (largest absolute eigenvalue)
    # of the shifted matrix.
    beta2, fields = linalg.eigs(A, k=order)
    
    # Find the mode with the lowest beta**2 out of those for which we solved.
    ind = np.argmin(beta2)

    # Return a WgMode class.
    return WgMode(np.sqrt(beta2[ind] - 4.0), fields[:, ind].real)
