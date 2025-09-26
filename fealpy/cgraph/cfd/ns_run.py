
from ..nodetype import CNodeType, PortConf, DataType

def lagrange_multiplier(A, b, c=0, uspace=None, pspace=None):
    """
    Constructs the augmented system matrix for Lagrange multipliers.
    c is the integral of pressure, default is 0.
    """
    from fealpy.sparse import COOTensor
    from fealpy.backend import backend_manager as bm
    from fealpy.fem import SourceIntegrator, LinearForm
    from fealpy.fem import BlockForm
    LagLinearForm = LinearForm(pspace)
    LagLinearForm.add_integrator(SourceIntegrator(source=1))
    LagA = LagLinearForm.assembly()
    LagA = bm.concatenate([bm.zeros(uspace.number_of_global_dofs()), LagA], axis=0)

    A1 = COOTensor(bm.array([bm.zeros(len(LagA), dtype=bm.int32),
                            bm.arange(len(LagA), dtype=bm.int32)]), LagA, spshape=(1, len(LagA)))

    A = BlockForm([[A, A1.T], [A1, None]])
    A = A.assembly_sparse_matrix(format='csr')
    b0 = bm.array([c])
    b  = bm.concatenate([b, b0], axis=0)
    return A, b

__all__ = ["StationaryNSRun"]

class StationaryNSRun(CNodeType):
    TITLE: str = "稳态 NS 方程有限元迭代求解"
    PATH: str = "流体.NS 方程有限元迭代求解"
    INPUT_SLOTS = [
        PortConf("maxstep", DataType.INT, 0, default=1000, min_val=1, title="最大迭代步数"),
        PortConf("tol", DataType.FLOAT, 0, default=1e-6, min_val=1e-12, max_val=1e-2, title="残差"),
        PortConf("update", DataType.FUNCTION, title="更新函数"),
        PortConf("apply_bc", DataType.FUNCTION, title="边界处理函数"),
        PortConf("BForm", DataType.LINOPS, title="算子"),
        PortConf("LForm", DataType.LINOPS, title="向量"),
        PortConf("uspace", DataType.SPACE, title="速度函数空间"),
        PortConf("pspace", DataType.SPACE, title="压力函数空间"),
        PortConf("mesh", DataType.MESH, title="网格")
    ]
    OUTPUT_SLOTS = [
        PortConf("uh", DataType.FUNCTION, title="速度数值解"),
        PortConf("uh_x", DataType.FUNCTION, title="速度x分量数值解"),
        PortConf("uh_y", DataType.FUNCTION, title="速度y分量数值解"),
        PortConf("ph", DataType.FUNCTION, title="压力数值解")
    ]
    @staticmethod
    def run(maxstep, tol, update, apply_bc, BForm, LForm, uspace, pspace, mesh):
        from fealpy.solver import spsolve
        uh0 = uspace.function()
        ph0 = pspace.function()
        uh1 = uspace.function()
        ph1 = pspace.function()
        ugdof = uspace.number_of_global_dofs()
        for i in range(maxstep):
            update(uh0)
            A = BForm.assembly()
            F = LForm.assembly()
            A, F = apply_bc(A, F)
            A, F = lagrange_multiplier(A, F, c = 0, uspace=uspace, pspace=pspace)
            x = spsolve(A, F,"mumps")
            uh1[:] = x[:ugdof]
            ph1[:] = x[ugdof:-1]
            res_u = mesh.error(uh0, uh1)
            res_p = mesh.error(ph0, ph1)
            
            if res_u + res_p < tol:
                break
            uh0[:] = uh1
            ph0[:] = ph1

        NN = mesh.number_of_nodes()
        uh_x = uh1[:int(ugdof/2)]
        uh_x = uh_x[:NN]
        uh_y = uh1[int(ugdof/2):]
        uh_y = uh_y[:NN]

        return uh1, uh_x, uh_y, ph1