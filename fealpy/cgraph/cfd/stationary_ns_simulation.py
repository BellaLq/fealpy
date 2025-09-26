
from ..nodetype import CNodeType, PortConf, DataType

__all__ = ["StationaryNSNewton", "StationaryNSOssen", "StationaryNSStokes"]


class StationaryNSNewton(CNodeType):
    TITLE: str = "稳态 NS 方程 Newton 迭代格式"
    PATH: str = "流体.有限元算法"
    INPUT_SLOTS = [
        PortConf("constitutive", DataType.MENU, 0, title="本构方程", default=1, items=[i for i in range(1, 2)]),
        PortConf("mu", DataType.FLOAT, title="粘度系数"),
        PortConf("rho", DataType.FLOAT, title = "密度"),
        PortConf("source", DataType.FUNCTION, title="源"),
        PortConf("uspace", DataType.SPACE, title="速度函数空间"),
        PortConf("pspace", DataType.SPACE, title="压力函数空间"),
        PortConf("q", DataType.INT, 0, default = 3, min_val=3, title="积分精度")
    ]
    OUTPUT_SLOTS = [
        PortConf("BForm", DataType.LINOPS, title="算子"),
        PortConf("LForm", DataType.LINOPS, title="向量"),
        PortConf("update", DataType.FUNCTION, title="更新函数")
    ]
    @staticmethod
    def run(constitutive, mu, rho, source, uspace, pspace, q):
        from fealpy.backend import backend_manager as bm
        from fealpy.decorator import barycentric
        from fealpy.fem import LinearForm, BilinearForm, BlockForm, LinearBlockForm
        from fealpy.fem import (ScalarMassIntegrator,ScalarConvectionIntegrator, PressWorkIntegrator, ScalarDiffusionIntegrator,
                                SourceIntegrator, ViscousWorkIntegrator)
 
        A00 = BilinearForm(uspace)
        u_BM_netwon = ScalarMassIntegrator(q = q)
        u_BC = ScalarConvectionIntegrator(q = q)
        if constitutive == 1:
            u_BVW = ScalarDiffusionIntegrator(q=q)
        elif constitutive == 2:
            u_BVW = ViscousWorkIntegrator(q=q)
        
        A00.add_integrator(u_BM_netwon)
        A00.add_integrator(u_BC)
        A00.add_integrator(u_BVW)
        
        A01 = BilinearForm((pspace, uspace))
        u_BPW = PressWorkIntegrator(q = q)
        A01.add_integrator(u_BPW)
       
        A = BlockForm([[A00, A01], [A01.T, None]]) 

        L0 = LinearForm(uspace)
        u_LSI = SourceIntegrator(q = q)
        u_source_LSI = SourceIntegrator(q = q)
        L0.add_integrator(u_LSI) 
        L0.add_integrator(u_source_LSI)
        L1 = LinearForm(pspace)
        L = LinearBlockForm([L0, L1])
        
        def update(u0): 
            cv = mu
            cc = rho
            pc = 1.0
            cbf = source
            
            ## BilinearForm
            u_BVW.coef = cv
            u_BPW.coef = -pc

            @barycentric
            def u_BC_coef(bcs, index):
                cccoef = cc(bcs, index)[..., bm.newaxis] if callable(cc) else cc
                return cccoef * u0(bcs, index)
            u_BC.coef = u_BC_coef

            @barycentric
            def u_BM_netwon_coef(bcs,index):
                cccoef = cc(bcs, index)[..., bm.newaxis] if callable(cc) else cc
                return cccoef * u0.grad_value(bcs, index)
            u_BM_netwon.coef = u_BM_netwon_coef

            ## LinearForm 
            @barycentric
            def u_LSI_coef(bcs, index):
                cccoef = cc(bcs, index)[..., bm.newaxis] if callable(cc) else cc
                result = cccoef*bm.einsum('...j, ...ij -> ...i', u0(bcs, index), u0.grad_value(bcs, index))
                return result
            u_LSI.source = u_LSI_coef
            u_source_LSI.source = cbf

        
        return (A, L, update)


class StationaryNSOssen(CNodeType):
    TITLE: str = "稳态 NS 方程 Ossen 迭代格式"
    PATH: str = "流体.有限元算法"
    INPUT_SLOTS = [
        PortConf("constitutive", DataType.MENU, 0, title="本构方程", default=1, items=[i for i in range(1, 2)]),
        PortConf("mu", DataType.FLOAT, title="粘度系数"),
        PortConf("rho", DataType.FLOAT, title = "密度"),
        PortConf("source", DataType.FUNCTION, title="源"),
        PortConf("uspace", DataType.SPACE, title="速度函数空间"),
        PortConf("pspace", DataType.SPACE, title="压力函数空间"),
        PortConf("q", DataType.INT, 0, default = 3, min_val=3, title="积分精度")
    ]
    OUTPUT_SLOTS = [
        PortConf("BForm", DataType.LINOPS, title="算子"),
        PortConf("LForm", DataType.LINOPS, title="向量"),
        PortConf("update", DataType.FUNCTION, title="更新函数")
    ]
    @staticmethod
    def run(constitutive, mu, rho, source, uspace, pspace, q):
        from fealpy.backend import backend_manager as bm
        from fealpy.decorator import barycentric
        from fealpy.fem import LinearForm, BilinearForm, BlockForm, LinearBlockForm
        from fealpy.fem import (ScalarConvectionIntegrator, PressWorkIntegrator, ScalarDiffusionIntegrator,
                                SourceIntegrator, ViscousWorkIntegrator)
 
        A00 = BilinearForm(uspace)
        u_BC = ScalarConvectionIntegrator(q=q)
        if constitutive == 1:
            u_BVW = ScalarDiffusionIntegrator(q=q)
        elif constitutive == 2:
            u_BVW = ViscousWorkIntegrator(q=q)
        
        A00.add_integrator(u_BC)
        A00.add_integrator(u_BVW)

        A01 = BilinearForm((pspace, uspace))
        u_BPW = PressWorkIntegrator(q=q)
        A01.add_integrator(u_BPW)
        
        A = BlockForm([[A00, A01], [A01.T, None]]) 

        L0 = LinearForm(uspace)
        u_source_LSI = SourceIntegrator(q=q)
        L0.add_integrator(u_source_LSI) 
        L1 = LinearForm(pspace)
        L = LinearBlockForm([L0, L1])

        def update(u0): 
            cv = mu
            cc = rho
            pc = 1.0
            cbf = source
            
            ## BilinearForm
            u_BVW.coef = cv
            u_BPW.coef = -pc

            @barycentric
            def u_BC_coef(bcs, index):
                cccoef = cc(bcs, index)[..., bm.newaxis] if callable(cc) else cc
                return cccoef * u0(bcs, index)
            u_BC.coef = u_BC_coef

            ## LinearForm 
            u_source_LSI.source = cbf

        return A, L, update



class StationaryNSStokes(CNodeType):
    TITLE: str = "稳态 NS 方程 Stokes 迭代格式"
    PATH: str = "流体.有限元算法"
    INPUT_SLOTS = [
        PortConf("constitutive", DataType.MENU, 0, title="本构方程", default=1, items=[i for i in range(1, 2)]),
        PortConf("mu", DataType.FLOAT, title="粘度系数"),
        PortConf("rho", DataType.FLOAT, title = "密度"),
        PortConf("source", DataType.FUNCTION, title="源"),
        PortConf("uspace", DataType.SPACE, title="速度函数空间"),
        PortConf("pspace", DataType.SPACE, title="压力函数空间"),   
        PortConf("q", DataType.INT, 0, default = 3, min_val=3, title="积分精度")
    ]
    OUTPUT_SLOTS = [
        PortConf("BForm", DataType.LINOPS, title="算子"),
        PortConf("LForm", DataType.LINOPS, title="向量"),
        PortConf("update", DataType.FUNCTION, title="更新函数")
    ]
    @staticmethod
    def run(constitutive, mu, rho, source, uspace, pspace, q):
        from fealpy.backend import backend_manager as bm
        from fealpy.fem import LinearForm, BilinearForm, BlockForm, LinearBlockForm
        from fealpy.fem import (PressWorkIntegrator, ScalarDiffusionIntegrator,
                            SourceIntegrator, ViscousWorkIntegrator)
        from fealpy.decorator import barycentric

        A00 = BilinearForm(uspace)
        
        if constitutive == 1:
            u_BVW = ScalarDiffusionIntegrator(q=q)
        elif constitutive == 2:
            u_BVW = ViscousWorkIntegrator(q=q)
        
        A00.add_integrator(u_BVW)
        
        A01 = BilinearForm((pspace, uspace))
        u_BPW = PressWorkIntegrator(q=q)
        A01.add_integrator(u_BPW)

        A = BlockForm([[A00, A01], [A01.T, None]]) 

        L0 = LinearForm(uspace)
        u_LSI = SourceIntegrator(q=q)
        u_source_LSI = SourceIntegrator(q=q)
        L0.add_integrator(u_LSI) 
        L0.add_integrator(u_source_LSI)
        L1 = LinearForm(pspace)
        L = LinearBlockForm([L0, L1])

        def update(u0): 
            cv = mu
            cc = rho
            pc = 1.0
            cbf = source
            
            ## BilinearForm
            u_BVW.coef = cv
            u_BPW.coef = -pc

            ## LinearForm 
            @barycentric
            def u_LSI_coef(bcs, index):
                cccoef = cc(bcs, index)[..., bm.newaxis] if callable(cc) else cc
                result = -cccoef*bm.einsum('...j, ...ij -> ...i', u0(bcs, index), u0.grad_value(bcs, index))
                return result
            u_LSI.source = u_LSI_coef
            u_source_LSI.source = cbf

        return A, L, update