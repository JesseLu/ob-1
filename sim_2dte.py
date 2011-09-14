import petsc4py
petsc4py.init('-ksp_monitor')
from petsc4py import PETSc
import numpy as np

class simulate:
    """ Used to solve Maxwell's equations in two-dimensions (TE mode).

    """

    def __init__(self, dims, omega, d_pml=5):
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

        stretch_coords(d_pml, dims[0]-0.5, pos, 1/omega)

        self.Dx_H = difference_matrix(self.da, (1, 0))
        self.Dy_H = difference_matrix(self.da, (0, 1))
        self.Dx_E = difference_matrix(self.da, (-1, 0))
        self.Dy_E = difference_matrix(self.da, (0, -1))

        # print self.Dx_H.getValues(range(9), range(9))
            
    def insert_structure(self):
        A = self.da.getMatrix() 
        A = self.Dx_E.matMult(self.Dx_H)
        A.axpy(1.0, self.Dy_E.matMult(self.Dy_H))
        A.shift(-self.omega**2)

        self.ksp = PETSc.KSP().create()
        self.ksp.setOperators(A)

    def solve(self, type, max_iters=1000):
        self.ksp.setFromOptions()
        self.ksp.setType(type)
        self.ksp.setTolerances(max_it=max_iters)
        # self.ksp.setMonitor(self.ksp.getMonitor())
        print self.ksp.getMonitor()
        self.ksp.solve(self.b, self.x)


def difference_matrix(da, shift, coeff=1.0, sc=None):
    """ Create a difference matrix.

    """

    D = da.getMatrix() # Create the difference matrix.

    row = PETSc.Mat.Stencil() # Defines row location of grid element.
    col = PETSc.Mat.Stencil() # Defines column location of grid element.

    (i0, i1), (j0, j1) = da.getRanges()
    for i in range(i0, i1):
        for j in range(j0, j1):
            row.index = (i, j)

            col.index = (i, j)
            D.setValueStencil(row, col, +1 * coeff)

            col.index = (i + shift[0], j + shift[1])
            D.setValueStencil(row, col, -1 * coeff)

    D.assemblyBegin()
    D.assemblyEnd()
    return D


def stretch_coords(d_pml, max_pos, pos, sigma, m=2.5):
    """ Generate the stretched coordinate coefficients which are used in
        the PML (perfectly matched layer) absorbing boundary.

    """

    sc = (1 + 0j) * np.ones(pos.shape)
    for k in range(pos.size):
        if pos[k] < d_pml:
            sc[k] = (1 + 1j * ((d_pml - pos[k]) / d_pml)**m)**-1
        elif max_pos - pos[k] < d_pml:
            sc[k] = (1 + 1j * ((d_pml - (max_pos - pos[k])) / d_pml)**m)**-1

    sc = (1 + 0j) * np.ones(pos.shape)
    return sc
