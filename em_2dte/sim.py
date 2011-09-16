import petsc4py
petsc4py.init('-ksp_type preonly -pc_type lu -ksp_monitor')
from petsc4py import PETSc
import numpy as np

class FieldSolve:
    """ Used to solve Maxwell's equations in two-dimensions (TE mode).

    """

    def __init__(self, dims, omega, d_pml=10):
        """ Create the matrices that we will use to simulate different 
            dielectric structures.

        """

        self.da = PETSc.DA().create(dims,
                                    stencil_width=1,
                                    boundary_type=('periodic', 'periodic'))

        self.omega = omega
        self.dims = dims
        self.b = self.da.createGlobalVec()
        self.x = self.da.createGlobalVec()

        pos = np.arange(dims[0])

        sc_x = lambda x, y: \
                stretch_coords(d_pml, dims[0]-0.5, x, 1/omega)
        sc_y = lambda x, y: \
                stretch_coords(d_pml, dims[1]-0.5, y, 1/omega)

        self.Dx_H = difference_matrix(self.da, (1, 0), 
                                        sc=lambda x, y: sc_x(x+0.5, y))
        self.Dy_H = difference_matrix(self.da, (0, 1),
                                        sc=lambda x, y: sc_y(x, y+0.5))
        self.Dx_E = difference_matrix(self.da, (-1, 0),
                                        sc=lambda x, y: sc_x(x, y))
        self.Dy_E = difference_matrix(self.da, (0, -1),
                                        sc=lambda x, y: sc_y(x, y))

        # print self.Dx_H.getValues(range(9), range(9))
            
    def insert_structure(self):
        A = self.da.getMatrix() 
        A = self.Dx_E.matMult(self.Dx_H)
        A.axpy(1.0, self.Dy_E.matMult(self.Dy_H))
        A.shift(-self.omega**2)

        self.ksp = PETSc.KSP().create()
        self.ksp.setOperators(A)

    def solve(self, max_iters=1000):
        self.ksp.setFromOptions()
#         self.ksp.setType('preonly')
#         print self.ksp.getPC()
#         self.ksp.setPC('lu')
        self.ksp.setTolerances(max_it=max_iters)
        print self.ksp.getType()
        # self.ksp.setMonitor(self.ksp.getMonitor())
        self.ksp.solve(self.b, self.x)


def difference_matrix(da, shift, coeff=1.0, sc=None):
    """ Create a difference matrix.

    """

    if sc is None:
        sc = lambda x, y: (1.0)

    D = da.getMatrix() # Create the difference matrix.

    row = PETSc.Mat.Stencil() # Defines row location of grid element.
    col = PETSc.Mat.Stencil() # Defines column location of grid element.

    (i0, i1), (j0, j1) = da.getRanges()
    for i in range(i0, i1):
        for j in range(j0, j1):
            row.index = (i, j)

            col.index = (i, j)
            D.setValueStencil(row, col, +1.0 * coeff * sc(i, j))

            col.index = (i + shift[0], j + shift[1])
            D.setValueStencil(row, col, -1.0 * coeff * sc(i, j))

    D.assemblyBegin()
    D.assemblyEnd()
    return D


def stretch_coords(d_pml, max_pos, pos, sigma, m=2.5):
    """ Generate the stretched coordinate coefficients which are used in
        the PML (perfectly matched layer) absorbing boundary.

    """

    if pos < d_pml:
        return (1 + 1j * ((d_pml - pos) / d_pml)**m)**-1
    elif max_pos - pos < d_pml:
        return (1 + 1j * ((d_pml - (max_pos - pos)) / d_pml)**m)**-1
    else:
        return 1

#     sc = (1 + 0j) * np.ones(pos.shape)
#     for k in range(pos.size):
#         if pos[k] < d_pml:
#             sc[k] = (1 + 1j * ((d_pml - pos[k]) / d_pml)**m)**-1
#         elif max_pos - pos[k] < d_pml:
#             sc[k] = (1 + 1j * ((d_pml - (max_pos - pos[k])) / d_pml)**m)**-1
# 
#     sc = (1 + 0j) * np.ones(pos.shape)
#     return sc
