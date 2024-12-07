from fealpy.backend import backend_manager as bm
from fealpy.functionspace import LagrangeFESpace,TensorFunctionSpace
from fealpy.fem import ScalarDiffusionIntegrator, BilinearForm, ScalarMassIntegrator, LinearForm
from fealpy.fem import ScalarSourceIntegrator, LinearElasticIntegrator,VectorSourceIntegrator
from fealpy.fem.dirichlet_bc import DirichletBC
from fealpy.solver import cg
from app.FuelRodSim.fuel_rod_mesher import FuelRodMesher
from fealpy.material.elastic_material import LinearElasticMaterial
from heat_elasticity_pde import Parabolic2dData
from heat_elasticity_pde import BoxDomainTriData2D

def varepsilon_therm(p0,p1,alpha_therm):
    delta_T = p1 - p0
    varepsilon_therm = alpha_therm * delta_T
    return varepsilon_therm

def sigma_eff(uh,material):
    uh = bm.array(uh)
    D = material.elastic_matrix()
    B = material.matrix()
    return sigma_eff
    
def varepsilon_sc(C,sigma_eff,phi,t):
    varepsilon_sc = C * sigma_eff * phi * (t * 3600)
    return varepsilon_sc

def varepsilon_irr():
    pass

#bm.set_backend('pytorch') # 选择后端为pytorch
mm = 1e-03
#包壳厚度
w = 0.15 * mm
#半圆半径
R1 = 0.5 * mm
#四分之一圆半径
R2 = 1.0 * mm
#连接处直线段
L = 0.575 * mm
#内部单元大小
h = 0.5 * mm
#棒长
l = 20 * mm

# 网格生成
mesher = FuelRodMesher(R1,R2,L,w,h,meshtype='segmented',modeltype='2D')
mesh = mesher.get_mesh
ficdx,cacidx = mesher.get_2D_fcidx_cacidx()
cnidx,bdnidx = mesher.get_2D_cnidx_bdnidx()

# 热传导
pde1=Parabolic2dData('exp(-2*pi**2*t)*sin(pi*x)*sin(pi*y)','x','y','t')
node = mesh.node
isBdNode = mesh.boundary_node_flag()
p0 = pde1.init_solution(node) #准备一个初值
#p = bm.array(p0)

space = LagrangeFESpace(mesh, p=1)
GD = space.geo_dimension()
duration = pde1.duration()
nt = 640
tau = (duration[1] -duration[0]) / nt


heat_alpha = 0.5
bform = BilinearForm(space)
bform.add_integrator(ScalarDiffusionIntegrator(heat_alpha, q=3))
heat_K = bform.assembly()

bform2 = BilinearForm(space)
bform2.add_integrator(ScalarMassIntegrator(q=3))
heat_M = bform2.assembly()
### 线弹性
pde2 = BoxDomainTriData2D()
tensor_space = TensorFunctionSpace(space, shape=(-1, 2))
gdof = space.number_of_global_dofs()
pfcm = LinearElasticMaterial(name='E1nu0.3', 
                                                elastic_modulus=1, poisson_ratio=0.3, 
                                                hypo='plane_strain', device=bm.get_device(mesh))
linear_integrator_K = LinearElasticIntegrator(material=pfcm, q=tensor_space.p+3)
bform = BilinearForm(tensor_space)
bform.add_integrator(linear_integrator_K)
linear_elasticity_K = bform.assembly(format='csr')

for n in range(nt):
    ### 热传导计算
    t = duration[0] + n * tau
    # 由于PDE模型基于符号计算，需要定义一个在笛卡尔坐标下的函数
    bform3 = LinearForm(space)
    from fealpy.decorator import cartesian
    @cartesian
    def coef(p):
        time = t
        val = pde1.source(p, time)
        return val
    bform3.add_integrator(ScalarSourceIntegrator(coef))
    F = bform3.assembly()
    A = heat_M +  heat_K * tau
    b = heat_M @ p0 + tau * F
    bc = DirichletBC(space=space,  gd=lambda p: pde1.dirichlet(p,t))
    A, b = bc.apply(A, b)
    p = cg(A, b, maxiter=5000, atol=1e-14, rtol=1e-14)
    ### 线弹性计算 #TODO 未完成 需要添加其他应变影响的项
    linear_integrator_F = VectorSourceIntegrator(source=pde2.source, q=tensor_space.p+3)
    #varepsilon = varepsilon_therm(p0,p,alpha)
    lform = LinearForm(tensor_space)    
    lform.add_integrator(linear_integrator_F)
    linear_elasticity_F = lform.assembly()
    dbc = DirichletBC(space=tensor_space, 
                    gd=pde2.dirichlet, 
                    threshold=None, 
                    method='interp')
    linear_elasticity_K, linear_elasticity_F = dbc.apply(A=linear_elasticity_K, f=linear_elasticity_F, uh=None, gd=pde2.dirichlet, check=True)
    uh = tensor_space.function()
    uh[:] = cg(linear_elasticity_K, linear_elasticity_F, maxiter=1000, atol=1e-14, rtol=1e-14)
    # 计算等效应力
    qf = mesh.quadrature_formula(q=tensor_space.p+3)
    bcs, ws = qf.get_quadrature_points_and_weights()
    gphi = tensor_space.grad_basis(bcs)
    NC = mesh.number_of_cells()
    tldof = tensor_space.number_of_local_dofs()
    cell2tdof = tensor_space.cell_to_dof()
    uh = bm.array(uh)  
    uh_cell = bm.zeros((NC, tldof)) # (NC, tldof)
    for c in range(NC):
        uh_cell[c] = uh[cell2tdof[c]]
    D = pfcm.elastic_matrix()
    B = pfcm.strain_matrix(True, gphi)
    # 重塑 u_h，使其适配 B 的形状
    uh_cell_new = uh_cell[:, None, :, None]  
    sigma = bm.einsum('ijkl, cmqijr, cmr -> cmqi', D, B, uh_cell_new)
    print(sigma.shape)

    p0 = p
        
    
    
    

