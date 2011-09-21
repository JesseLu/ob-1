import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg

class WgMode:
    """ Class for working with waveguide modes.

    """

    def __init__(self, beta, field, error):

        self.beta = beta # Wave-vector of the mode.
        self.field = field # Field profile of the mode.
        self.error = error # Numerical error in eigenvalue equation.

        # If beta is real (i.e. we have a true propagating mode), 
        # then make sure that the fields are all real as well.
        if np.isreal(self.beta):
            if all(np.isreal(self.field)):
                self.beta = np.real(self.beta)
                self.field = np.real(self.field)
            else:
                raise Exception("Waveguide mode should be real-valued.")

        # Orient the field so that the real part of it's largest component is
        # positive.
        ind = np.argmax(np.abs(self.field))
        if np.real(self.field[ind]) < 0:
            self.field = -1 * self.field

        # Note that many modes will have symmetry such that their positive
        # and negative real parts are equal, in this case, an additional
        # standard is needed.

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
    epsilon_x = 0.5 * (epsilon + epsilon[range(1, n) + [0]])
    # epsilon_x = 0.5 * (epsilon + epsilon[[n-1] + range(n-1)])
    epsilon_y = epsilon

    # Helper function to quickly create a sparse, diagonal matrix.
    my_diag = lambda x: sparse.spdiags(x, 0, n, n)

    # Form A.
    A = -1 * my_diag(epsilon_x) * Dx.T * my_diag(epsilon_y**-1) * Dx + \
        omega**2 * my_diag(epsilon_x)

    # Calculate relevant eigenvalues.
    # *   These will be the maximal eigenvalues (largest absolute eigenvalue)
    #     of the shifted matrix. Shifting value of 4 is used since most 
    #     negative beta possible is -2 (and 4 = (-2)^2).
    # *   Force atleast 5 eigenvalues to be found, for some reason, the 
    #     underlying linear algebra package doesn't like to find only 1 or 2 
    #     eigenvalues.
    beta2, fields = linalg.eigs(A + 4*sparse.eye(n, n), k=np.max([order, 5]))
        
    # Sort modes by eigenvalue. Used to extract the mode of interest.
    ind = np.argsort(beta2)

    # Pick out the mode of interest
    beta = np.sqrt(beta2[ind[-order]] - 4.0)
    field = fields[:, ind[-order]]

    # Calculate the error
    error = np.linalg.norm(A * field - beta**2 * field)

    # Return a WgMode class containing the correct waveguide mode.
    return WgMode(beta, field, error)
