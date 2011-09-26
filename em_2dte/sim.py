import petsc4py
from petsc4py import PETSc
import numpy as np
import wgmode

class FieldSolve:
    """ Holds the result of a simulation of Maxwell's equations.

    """

    def __init__(self, omega, mode_order, epsilon, field, res_norm):
        """ Store the relevant parameters of the simulation.

        """
        self.omega = omega
        self.mode_order = mode_order
        self.epsilon = epsilon
        self.field = field 
        self.res_norm = res_norm


    def plot(self, gplot):
        gplot.plot(np.real(self.field))
        # gplot.plot(np.abs(self.field))

def solve(omega, epsilon, mode_pos=10, mode_order=1, d_pml=10.):

    # Create the distributed array that all vectors and matrices will be 
    # based off of.
    da = PETSc.DA().create( epsilon.shape,
                            stencil_width=1,
                            boundary_type=('periodic', 'periodic'))

    # Create the matrix used to solve for the field.
    A = _setup_matrix(da, epsilon, omega, d_pml)

    # Determine current sources needed to excite a 1-way waveguide mode.
    mode = wgmode.solve(omega, epsilon[mode_pos,:], mode_order) # Find mode.

    b = da.createGlobalVec() # Create current source vector.
    b_val = da.getVecArray(b) # Obtain access to elements of b.
    b_val[mode_pos, :] = mode.field[:] # Insert mode at two locations.
    b_val[mode_pos-1, :] = -mode.field[:] * np.exp(1j * mode.beta) 


    # Solve the linear system.
    ksp = PETSc.KSP().create()
    ksp.setOperators(A)


    ksp.setType('preonly') # Don't use any Krylov subspace iterative algorithms.
    pc = ksp.getPC()
    pc.setType('lu') # Instead, solve using good ol' LU factorization.

#     ksp.setType('bcgs') # Use Bi-CG-Stab algorithm.
#     pc = ksp.getPC()
#     pc.setType('none') # Without preconditioning.
#     ksp.setTolerances(max_it=10000) # Set maximum iterations to 10,000.

    # Note that I should add monitoring to this, but haven't been able to 
    # figure it out (2011-09-26).

    x = da.createGlobalVec() # Vector to store solution in.
    ksp.solve(b, x) # Solve.

    # Print result.
    print "Simulated in", ksp.getIterationNumber(), "iterations,", \
            "with residual norm of", ksp.getResidualNorm()

    # Return result as a FieldSolve class instance.
    return FieldSolve(omega, mode_order, epsilon, 
                        da.getVecArray(x)[:], ksp.getResidualNorm())


def _setup_matrix(da, epsilon, omega, d_pml):
    m = np.sqrt(np.max(epsilon.flatten()))

    sc_x = lambda x, y: \
            _stretch_coords(d_pml, epsilon.shape[0]-0.5, x, 1/(omega * m))
    sc_y = lambda x, y: \
            _stretch_coords(d_pml, epsilon.shape[1]-0.5, y, 1/(omega * m))

    # H-fields need to "reach" back since that's where corresponding E-fields
    # are located.
    Dx_H = _diff_mat(da, (-1, 0), stretch_coords=lambda x, y: sc_x(x, y))
    Dy_H = _diff_mat(da, (0, -1), stretch_coords=lambda x, y: sc_y(x, y))

    Dx_E = _diff_mat(da, (1, 0), stretch_coords=lambda x, y: sc_x(x+0.5, y))
    Dy_E = _diff_mat(da, (0, 1), stretch_coords=lambda x, y: sc_y(x, y+0.5))

    inv_ex, inv_ey = _get_inv_eps(da, epsilon)

    Dx_H.diagonalScale(inv_ey)
    A = Dx_E.matMult(Dx_H)

    Dy_H.diagonalScale(inv_ex)
    A.axpy(1.0, Dy_E.matMult(Dy_H))

    A.shift(-omega**2)

    return A


def _diff_mat(da, shift, coeff=1.0, stretch_coords=None):
    """ Create a difference matrix.

    """

    if stretch_coords is None:
        stretch_coords = lambda x, y: (1.0)

    D = da.getMatrix() # Create the difference matrix.

    row = PETSc.Mat.Stencil() # Defines row location of grid element.
    col = PETSc.Mat.Stencil() # Defines column location of grid element.

    (i0, i1), (j0, j1) = da.getRanges()
    for i in range(i0, i1):
        for j in range(j0, j1):
            row.index = (i, j)

            col.index = (i, j)
            D.setValueStencil(row, col, +1.0 * coeff / stretch_coords(i, j))

            col.index = (i + shift[0], j + shift[1])
            D.setValueStencil(row, col, -1.0 * coeff / stretch_coords(i, j))

    D.assemblyBegin()
    D.assemblyEnd()
    return D


def _stretch_coords(d_pml, max_pos, pos, sigma, m=2.5):
    """ Generate the stretched coordinate coefficients which are used in
        the PML (perfectly matched layer) absorbing boundary.

    """

    if pos < d_pml:
        return (1 + 1j * sigma * (float(d_pml - pos) / d_pml)**m)
    elif max_pos - pos < d_pml:
        return (1 + 1j * sigma * (float(d_pml - (max_pos - pos)) / d_pml)**m)
    else:
        return 1

def _get_inv_eps(da, epsilon):
    shape = epsilon.shape

    inv_ex = da.createGlobalVec()
    da.getVecArray(inv_ex)[:] = (0.5 * (epsilon[:,:] + 
                            epsilon[range(1, shape[0]) + [0], :]))**-1

    inv_ey = da.createGlobalVec()
    da.getVecArray(inv_ey)[:] = (0.5 * (epsilon[:,:] +
                            epsilon[:, range(1, shape[1]) + [0]]))**-1

    return inv_ex, inv_ey
