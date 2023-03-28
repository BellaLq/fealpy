import numpy as np


class TrussStructureIntegrator:
    def __init__(self, E, A):
        self.E = E  # 杨氏模量
        self.A = A  # 单元横截面积

    def assembly_cell_matrix(self, space, index=np.s_[:], cellmeasure=None):
        """
        @brief 
        """
        assert isinstance(space, tuple) 
        space0 = space[0]
        mesh = space0.mesh
        GD = mesh.geo_dimension()

        assert  len(space) == GD

        c = self.E*self.A
        NC = mesh.number_of_cells()
        l = mesh.entity_measure('cell')
        tan = mesh.cell_unit_tangent()

        R = np.einsum('ik, im, i->ikm', tan, tan)
        R *=c/l[:, None, None]

        K = np.zeros((NC, 2*GD, 2*GD), dtype=np.float64)
        K[:, :GD, :GD] = R
        K[:, -GD:, :GD] = -R
        K[:, :GD, -GD:] = -R
        K[:, -GD:, -GD:] = R

        return (R, -R), (-R, R) 