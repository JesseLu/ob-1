import petsc4py
petsc4py.init('-ksp_type preonly -pc_type lu -ksp_monitor')
from petsc4py import PETSc
import numpy as np
import wgmode

class FieldSolve:
    """ Used to solve Maxwell's equations in two-dimensions (TE mode).

    """

    def __init__(self, field):
        """ Create the matrices that we will use to simulate different 
            dielectric structures.

        """
        self.field = field 



    def solve(self, omega, max_iters=1000):
        self.A.shift(-omega**2)

        self.ksp = PETSc.KSP().create()
        self.ksp.setOperators(self.A)
        self.ksp.setFromOptions()
#         self.ksp.setType('preonly')
#         print self.ksp.getPC()
#         self.ksp.setPC('lu')
        self.ksp.setTolerances(max_it=max_iters)
        print self.ksp.getType()
        # self.ksp.setMonitor(self.ksp.getMonitor())
        self.x = self.da.createGlobalVec()
        self.ksp.solve(self.b, self.x)

    def plot(self, gplot):
        gplot.plot(np.real(self.field))
        gplot.plot(np.abs(self.field))

def solve(omega, epsilon, mode_pos=20, mode_order=1, d_pml=10):

    # Create the distributed array that all vectors and matrices will be 
    # based off of.
    da = PETSc.DA().create( epsilon.shape,
                            stencil_width=1,
                            boundary_type=('periodic', 'periodic'))

    # Create the matrix used to solve for the field.
    A = setup_matrix(da, epsilon, omega, d_pml)

    # Determine current sources
    mode = wgmode.solve(omega, epsilon[mode_pos,:], mode_order)
    b = da.createGlobalVec()
    b_val = da.getVecArray(b) # Obtain access to elements of b.
    b_val[mode_pos, :] = mode.field[:] # Define location of point source.
#     b_val[mode_pos+1, :] = mode.field[:] * np.exp(1j * np.pi) # Define location of point source.
#     b_val[mode_pos+1, :] = mode.field[:] * np.exp(1j * 2 * np.pi / 9) # Define location of point source.
    print mode.beta, 2*np.pi /mode.beta # Define location of point source.


    ksp = PETSc.KSP().create()
    ksp.setOperators(A)
    ksp.setFromOptions()
#         ksp.setType('preonly')
#         print ksp.getPC()
#         ksp.setPC('lu')
    ksp.setTolerances(max_it=100)
    print ksp.getType()
    # ksp.setMonitor(ksp.getMonitor())
    x = da.createGlobalVec()
    ksp.solve(b, x)

    return FieldSolve(da.getVecArray(x)[:])


def setup_matrix(da, epsilon, omega, d_pml):
    m = np.sqrt(np.max(epsilon.flatten()))

    sc_x = lambda x, y: \
            stretch_coords(d_pml, epsilon.shape[0]-0.5, x, 1/(omega * m))
    sc_y = lambda x, y: \
            stretch_coords(d_pml, epsilon.shape[1]-0.5, y, 1/(omega * m))

    # H-fields need to "reach" back since that's where corresponding E-fields
    # are located.
    Dx_H = diff_mat(da, (-1, 0), stretch_coords=lambda x, y: sc_x(x+0.5, y))
    Dy_H = diff_mat(da, (0, -1), stretch_coords=lambda x, y: sc_y(x, y+0.5))

    Dx_E = diff_mat(da, (1, 0), stretch_coords=lambda x, y: sc_x(x, y))
    Dy_E = diff_mat(da, (0, 1), stretch_coords=lambda x, y: sc_y(x, y))

    ex, ey = interp_epsilon(da, epsilon)

    Dx_H.diagonalScale(ey)
    A = Dx_E.matMult(Dx_H)

    Dy_H.diagonalScale(ex)
    A.axpy(1.0, Dy_E.matMult(Dy_H))

    A.shift(-omega**2)

    return A


def diff_mat(da, shift, coeff=1.0, stretch_coords=None):
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


def stretch_coords(d_pml, max_pos, pos, sigma, m=2.5):
    """ Generate the stretched coordinate coefficients which are used in
        the PML (perfectly matched layer) absorbing boundary.

    """

    if pos < d_pml:
        return (1 + 1j * sigma * ((d_pml - pos) / d_pml)**m)**1
    elif max_pos - pos < d_pml:
        return (1 + 1j * sigma * ((d_pml - (max_pos - pos)) / d_pml)**m)**1
    else:
        return 1

def interp_epsilon(da, epsilon):
    shape = epsilon.shape

    ex = da.createGlobalVec()
    da.getVecArray(ex)[:] = 0.5 * (epsilon[:,:] + 
                            epsilon[range(1, shape[0]) + [0], :])**-1

    ey = da.createGlobalVec()
    da.getVecArray(ey)[:] = 0.5 * (epsilon[:,:] +
                            epsilon[:, range(1, shape[1]) + [0]])**-1

    return ex, ey
