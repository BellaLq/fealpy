import taichi as ti
import numpy as np
from scipy.sparse import csr_matrix


def construct_edge(cell):
    """ 
    """
    NC =  cell.shape[0] 
    NEC = 3 
    NVE = 2 

    localEdge = np.array([(1, 2), (2, 0), (0, 1)])
    totalEdge = cell[:, localEdge].reshape(-1, 2)
    _, i0, j = np.unique(np.sort(totalEdge, axis=-1),
            return_index=True,
            return_inverse=True,
            axis=0)
    NE = i0.shape[0]
    edge2cell = np.zeros((NE, 4), dtype=np.int_)

    i1 = np.zeros(NE, dtype=np.int_)
    i1[j] = np.arange(NEC*NC, dtype=np.int_)

    edge2cell[:, 0] = i0//3
    edge2cell[:, 1] = i1//3
    edge2cell[:, 2] = i0%3
    edge2cell[:, 3] = i1%3

    edge = totalEdge[i0, :]
    cell2edge = np.reshape(j, (NC, 3))
    return edge, edge2cell, cell2edge


@ti.data_oriented
class TriangleMesh():
    def __init__(self, node, cell, itype=ti.u32, ftype=ti.f64):
        assert cell.shape[-1] == 3

        self.itype = itype
        self.ftype = ftype

        NN = node.shape[0]
        GD = node.shape[1]
        self.node = ti.field(self.ftype, (NN, GD))
        self.node.from_numpy(node)

        NC = cell.shape[0]
        self.cell = ti.field(self.itype, shape=(NC, 3))
        self.cell.from_numpy(cell)

        self.glambda = ti.field(self.ftype, shape=(NC, 3, GD)) 
        self.cellmeasure = ti.field(self.ftype, shape=(NC, ))
        self.init_grad_lambdas()

        self.construct_data_structure(cell)

    def construct_data_structure(self, cell):

        edge, edge2cell, cell2edge = construct_edge(cell)
        NE = edge.shape[0]

        self.edge = ti.field(self.itype, shape=(NE, 2))
        self.edge.from_numpy(edge)

        self.edge2cell = ti.field(self.itype, shape=(NE, 4))
        self.edge2cell.from_numpy(edge2cell)

        NC = self.number_of_cells()
        self.cell2edge = ti.field(self.itype, shape=(NC, 3))
        self.cell2edge.from_numpy(cell2edge)

    def multi_index_matrix(self, p):
        ldof = (p+1)*(p+2)//2
        idx = np.arange(0, ldof)
        idx0 = np.floor((-1 + np.sqrt(1 + 8*idx))/2)
        multiIndex = np.zeros((ldof, 3), dtype=np.int_)
        multiIndex[:,2] = idx - idx0*(idx0 + 1)/2
        multiIndex[:,1] = idx0 - multiIndex[:,2]
        multiIndex[:,0] = p - multiIndex[:, 1] - multiIndex[:, 2]
        return multiIndex

    def number_of_local_interpolation_points(self, p):
        return (p+1)*(p+2)//2

    def number_of_global_interpolation_points(self, p):
        NP = self.number_of_nodes()
        if p > 1:
            NE = self.number_of_edges()
            NP += (p-1)*NE
        if p > 2:
            NC = self.number_of_cells()
            NP += (p-2)*(p-1)*NC//2
        return NP

    @ti.kernel 
    def interpolation_points(self, p: ti.u32, ipoints: ti.template()):
        NN = self.node.shape[0]
        GD = self.node.shape[1]
        for i in range(NN):
            for j in range(GD):
                ipoints[i, j] = self.node[i, j]
        if p > 1:
            NE = self.edge.shape[0]
            for i in range(NE):
                for j in range(0, p-1):
                    I = NN + i*(p-1) + j
                    for k in range(GD):
                        ipoints[I, k] = ((p-j-1)*self.node[self.edge[i, 0], k] + (j+1)*self.node[self.edge[i, 0], k])/p
        if p > 2:
            NC = self.cell.shape[0]
            for i in range(NC):
                pass

    @ti.kernel
    def cell_to_dof(self, p: ti.u32, cell2dof: ti.template()):
        cdof = (p+1)*(p+2)//2 
        NN = self.node.shape[0]
        NE = self.edge.shape[0]
        for c in range(self.cell.shape[0]): 
            # 三个顶点 
            cell2dof[c, 0] = self.cell[c, 0]
            cell2dof[c, cdof - p - 1] = self.cell[c, 1] # 不支持负数索引
            cell2dof[c, cdof - 1] = self.cell[c, 2]

            # 第 0 条边
            e = self.cell2edge[c, 0]
            v0 = self.edge[e, 0]
            s0 = cdof - p
            if v0 == self.cell[c, 1]:
                for i in range(0, p-1):
                    cell2dof[c, s0] = NN + e*(p-1) + i 
                    s0 += 1
            else:
                for i in range(0, p-1):
                    cell2dof[c, s0] = NN + e*(p-1) + p - 2 - i
                    s0 += 1

            # 第 1 条边
            e = self.cell2edge[c, 1]
            v0 = self.edge[e, 0]
            s0 = 2
            if v0 == self.cell[c, 0]:
                for i in range(0, p-1):
                    cell2dof[c, s0] = NN + e*(p-1) + i 
                    s0 += i + 2
            else:
                for i in range(0, p-1):
                    cell2dof[c, s0] = NN + e*(p-1) + p - 2 - i 
                    s0 += i + 2 

            # 第 2 条边
            e = self.cell2edge[c, 2]
            v0 = self.edge[e, 0]
            s0 = 1
            if v0 == self.cell[c, 0]:
                for i in range(0, p-1):
                    cell2dof[c, s0] = NN + e*(p-1) + i 
                    s0 += i + 1
            else:
                for i in range(0, p-1):
                    cell2dof[c, s0] = NN + e*(p-1) + p - 2 - i 
                    s0 += i + 1 

            # 内部点
            if p >= 3:
                start = NN + (p-1)*NE
                level = p - 2 
                nd = (p - 2)*(p - 1)//2
                s0 = 4
                s1 = 0
                for l in range(0, level):
                    for i in range(0, l+1):
                        cell2dof[c, s0] = start + c*nd + s1 
                        s0 += 1
                        s1 += 1 
                    s0 += 3

    def number_of_nodes(self):
        return self.node.shape[0]

    def number_of_edges(self):
        return self.edge.shape[0]

    def number_of_cells(self):
        return self.cell.shape[0]

    def geo_dimension(self):
        return self.node.shape[1]

    def top_dimension(self):
        return 2 

    def entity(self, etype=2):
        if etype in {'cell', 2}:
            return self.cell
        elif etype in {'edge', 'face', 1}:
            return self.edge
        elif etype in {'node', 0}:
            return self.node
        else:
            raise ValueError("`etype` is wrong!")

    def vtk_cell_type(self, etype='cell'):
        if etype in {'cell', 2}:
            VTK_TRIANGLE = 5
            return VTK_TRIANGLE
        elif etype in {'face', 'edge', 1}:
            VTK_LINE = 3
            return VTK_LINE

    def to_vtk(self, fname, nodedata=None, celldata=None):
        """
        Parameters
        ----------

        Notes
        -----
        把网格转化为 VTK 的格式
        """
        from fealpy.mesh.vtk_extent import vtk_cell_index, write_to_vtu

        node = self.node.to_numpy()
        GD = self.geo_dimension()
        if GD == 2:
            node = np.concatenate((node, np.zeros((node.shape[0], 1))), axis=1)

        cell = self.cell.to_numpy(dtype=np.int_)
        cellType = self.vtk_cell_type()
        NV = cell.shape[-1]

        cell = np.r_['1', np.zeros((len(cell), 1), dtype=cell.dtype), cell]
        cell[:, 0] = NV

        NC = len(cell)
        print("Writting to vtk...")
        write_to_vtu(fname, node, NC, cellType, cell.flatten(),
                nodedata=nodedata,
                celldata=celldata)

    @ti.kernel
    def init_grad_lambdas(self):
        """
        初始化网格中每个单元上重心坐标函数的梯度，以及单元的面积
        """
        assert self.node.shape[1] == 2

        for i in range(self.cell.shape[0]):
            x0 = self.node[self.cell[i, 0], 0]
            y0 = self.node[self.cell[i, 0], 1]

            x1 = self.node[self.cell[i, 1], 0]
            y1 = self.node[self.cell[i, 1], 1]

            x2 = self.node[self.cell[i, 2], 0]
            y2 = self.node[self.cell[i, 2], 1]

            l = (x1 - x0)*(y2 - y0) - (y1 - y0)*(x2 - x0) 

            self.cellmeasure[i] = 0.5*l
            self.glambda[i, 0, 0] = (y1 - y2)/l
            self.glambda[i, 0, 1] = (x2 - x1)/l 
            self.glambda[i, 1, 0] = (y2 - y0)/l 
            self.glambda[i, 1, 1] = (x0 - x2)/l
            self.glambda[i, 2, 0] = (y0 - y1)/l
            self.glambda[i, 2, 1] = (x1 - x0)/l


    @ti.func
    def cell_measure(self, i: ti.u32) -> ti.f64:
        x0 = self.node[self.cell[i, 0], 0]
        y0 = self.node[self.cell[i, 0], 1]

        x1 = self.node[self.cell[i, 1], 0]
        y1 = self.node[self.cell[i, 1], 1]

        x2 = self.node[self.cell[i, 2], 0]
        y2 = self.node[self.cell[i, 2], 1]

        l = (x1 - x0)*(y2 - y0) - (y1 - y0)*(x2 - x0) 
        l *= 0.5
        return l

    @ti.func
    def grad_lambda(self, i: ti.u32) -> (ti.types.matrix(3, 2, ti.f64), ti.f64):
        """
        计算第 i 个单元上重心坐标函数的梯度，以及单元的面积
        """

        assert self.node.shape[1] == 2

        x0 = self.node[self.cell[i, 0], 0]
        y0 = self.node[self.cell[i, 0], 1]

        x1 = self.node[self.cell[i, 1], 0]
        y1 = self.node[self.cell[i, 1], 1]

        x2 = self.node[self.cell[i, 2], 0]
        y2 = self.node[self.cell[i, 2], 1]

        l = (x1 - x0)*(y2 - y0) - (y1 - y0)*(x2 - x0) 

        gphi = ti.Matrix.zero(ti.f64, 3, 2)
        gphi[0, 0] = (y1 - y2)/l
        gphi[0, 1] = (x2 - x1)/l 
        gphi[1, 0] = (y2 - y0)/l 
        gphi[1, 1] = (x0 - x2)/l
        gphi[2, 0] = (y0 - y1)/l
        gphi[2, 1] = (x1 - x0)/l

        l *= 0.5
        return gphi, l

    @ti.func
    def surface_grad_lambda(self, i: ti.u32) -> (ti.types.matrix(3, 2, ti.f64), ti.f64):
        """
        计算第 i 个单元上重心坐标函数的梯度，以及单元的面积
        """
        assert self.node.shape[1] == 3

        x0 = self.node[self.cell[i, 0], 0]
        y0 = self.node[self.cell[i, 0], 1]
        z0 = self.node[self.cell[i, 0], 2]

        x1 = self.node[self.cell[i, 1], 0]
        y1 = self.node[self.cell[i, 1], 1]
        z1 = self.node[self.cell[i, 0], 2]

        x2 = self.node[self.cell[i, 2], 0]
        y2 = self.node[self.cell[i, 2], 1]
        z2 = self.node[self.cell[i, 0], 2]

        gphi = ti.Matrix.zero(ti.f64, 3, 3)
        return grad, l

    @ti.kernel
    def cell_stiff_matrices(self, S: ti.template()):
        """
        计算网格上的所有单元刚度矩阵
        """
        for c in range(self.cell.shape[0]):
            gphi, cm = self.grad_lambda(c) 

            S[c, 0, 0] = cm*(gphi[0, 0]*gphi[0, 0] + gphi[0, 1]*gphi[0, 1])
            S[c, 0, 1] = cm*(gphi[0, 0]*gphi[1, 0] + gphi[0, 1]*gphi[1, 1])
            S[c, 0, 2] = cm*(gphi[0, 0]*gphi[2, 0] + gphi[0, 1]*gphi[2, 1])

            S[c, 1, 0] = S[c, 0, 1]
            S[c, 1, 1] = cm*(gphi[1, 0]*gphi[1, 0] + gphi[1, 1]*gphi[1, 1])
            S[c, 1, 2] = cm*(gphi[1, 0]*gphi[2, 0] + gphi[1, 1]*gphi[2, 1])

            S[c, 2, 0] = S[c, 0, 2]
            S[c, 2, 1] = S[c, 1, 2]
            S[c, 2, 2] = cm*(gphi[2, 0]*gphi[2, 0] + gphi[2, 1]*gphi[2, 1])

    @ti.kernel
    def cell_mass_matrices(self, S: ti.template()):
        """
        计算网格上的所有单元质量矩阵
        """
        for c in range(self.cell.shape[0]):
            cm = self.cell_measure(c)
            c0 = cm/6.0
            c1 = cm/12.0

            S[c, 0, 0] = c0 
            S[c, 0, 1] = c1
            S[c, 0, 2] = c1

            S[c, 1, 0] = c1 
            S[c, 1, 1] = c0  
            S[c, 1, 2] = c1

            S[c, 2, 0] = c1 
            S[c, 2, 1] = c1 
            S[c, 2, 2] = c0 

    @ti.kernel
    def cell_convection_matrices(self, u: ti.template(), S: ti.template()):
        """
        计算网格上所有单元的对流矩阵
        """
        for c in range(self.cell.shape[0]):
            gphi, cm = self.grad_lambda(c) 

            c0 = cm/6.0
            c1 = cm/12.0

            U = ti.Matrix.zero(ti.f64, 3, 2)
            for i in ti.static(range(2)):
                U[0, i] += u[self.cell[c, 0], i]*c0 
                U[0, i] += u[self.cell[c, 1], i]*c1 
                U[0, i] += u[self.cell[c, 2], i]*c1 

            for i in ti.static(range(2)):
                U[1, i] += u[self.cell[c, 0], i]*c1 
                U[1, i] += u[self.cell[c, 1], i]*c0 
                U[1, i] += u[self.cell[c, 2], i]*c1 

            for i in ti.static(range(2)):
                U[2, i] += u[self.cell[c, 0], i]*c1 
                U[2, i] += u[self.cell[c, 1], i]*c1 
                U[2, i] += u[self.cell[c, 2], i]*c0 

            for i in ti.static(range(3)):
                for j in ti.static(range(3)):
                    S[c, i, j] = U[i, 0]*gphi[j, 0] + U[i, 1]*gphi[j, 1]


    @ti.kernel
    def cell_source_vectors(self, f:ti.template(), bc:ti.template(), ws:ti.template(), F:ti.template()):
        """
        计算所有单元载荷
        """
        for c in range(self.cell.shape[0]):
            x0 = self.node[self.cell[c, 0], 0]
            y0 = self.node[self.cell[c, 0], 1]

            x1 = self.node[self.cell[c, 1], 0]
            y1 = self.node[self.cell[c, 1], 1]

            x2 = self.node[self.cell[c, 2], 0]
            y2 = self.node[self.cell[c, 2], 1]
            l = (x1 - x0)*(y2 - y0) - (y1 - y0)*(x2 - x0) 
            l *= 0.5
            for q in ti.static(range(bc.n)):
                x = x0*bc[q, 0] + x1*bc[q, 1] + x2*bc[q, 1]
                y = y0*bc[q, 0] + y1*bc[q, 1] + y2*bc[q, 1]
                z = f(x, y)
                for i in ti.static(range(3)):
                    F[c, i] += ws[q]*bc[q, i]*z

            for i in range(3):
                F[c, i] *= l

    @ti.kernel
    def ti_cell_stiff_matrices(self, K: ti.types.sparse_matrix_builder()):
        """
        
        """
        for c in range(self.cell.shape[0]):
            gphi, cm = self.grad_lambda(c) 
            for i in ti.static(range(3)):
                for j in ti.static(range(3)):
                    I = self.cell[c, i]
                    J = self.cell[c, j]
                    K[I, J] += cm*(gphi[i, 0]*gphi[j, 0] + gphi[i, 1]*gphi[j, 1]) 

    def ti_stiff_matrix(self, c=None):
        """
        基于 Taichi 组装刚度矩阵
        """
        NN = self.number_of_nodes()
        NC = self.number_of_cells()
        K = ti.linalg.SparseMatrixBuilder(NN, NN, max_num_triplets=9*NC)
        self.ti_cell_stiff_matrices(K)
        A = K.build()
        return A

    def stiff_matrix(self, c=None):
        """
        组装总体刚度矩阵
        """
        NC = self.number_of_cells()
        K = ti.field(ti.f64, (NC, 3, 3))
        self.cell_stiff_matrices(K)

        M = K.to_numpy()
        if c is not None:
            M *= c # 目前假设 c 为一常数

        cell = self.cell.to_numpy()
        I = np.broadcast_to(cell[:, :, None], shape=M.shape)
        J = np.broadcast_to(cell[:, None, :], shape=M.shape)

        NN = self.number_of_nodes()
        M = csr_matrix((K.to_numpy().flat, (I.flat, J.flat)), shape=(NN, NN))
        return M

    def mass_matrix(self, c=None):
        """
        组装总体质量矩阵
        """
        NC = self.number_of_cells() 

        K = ti.field(ti.f64, (NC, 3, 3))
        self.cell_mass_matrices(K)

        M = K.to_numpy()
        if c is not None:
            M *= c

        cell = self.cell.to_numpy()
        I = np.broadcast_to(cell[:, :, None], shape=M.shape)
        J = np.broadcast_to(cell[:, None, :], shape=M.shape)

        NN = self.number_of_nodes() 
        M = csr_matrix((M.flat, (I.flat, J.flat)), shape=(NN, NN))
        return M

    def convection_matrix(self, u):
        """
        组装总体对流矩阵
        """

        NC = self.number_of_cells()
        C = ti.field(ti.f64, (NC, 3, 3))
        self.cell_convection_matrices(u, C)

        M = C.to_numpy()

        cell = self.cell.to_numpy()
        I = np.broadcast_to(cell[:, :, None], shape=M.shape)
        J = np.broadcast_to(cell[:, None, :], shape=M.shape)

        NN = self.number_of_nodes()
        M = csr_matrix((M.flat, (I.flat, J.flat)), shape=(NN, NN))
        return M


    def source_vector(self, f):
        """
        组装总体载荷向量
        """
        NN = self.node.shape[0]
        NC = self.cell.shape[0]
        bc = ti.Matrix([
            [0.6666666666666670,	0.1666666666666670,     0.1666666666666670],
            [0.1666666666666670,	0.6666666666666670,     0.1666666666666670],
            [0.1666666666666670,	0.1666666666666670, 0.6666666666666670]], dt=ti.f64)
        ws = ti.Vector([0.3333333333333330, 0.3333333333333330, 0.3333333333333330], dt=ti.f64)

        F = ti.field(ti.f64, (NC, 3))
        self.cell_source_vectors(f, bc, ws, F)

        bb = F.to_numpy()
        F = np.zeros(NN, dtype=np.float64)
        cell = self.cell.to_numpy()
        np.add.at(F, cell, bb)
        return F

    @ti.kernel
    def scalar_linear_interpolation(self, f: ti.template(), R: ti.template()):
        """! 定义在二维平面上的标量函数的线性插值
        """
        for i in range(self.node.shape[0]):
            R[i] = f(self.node[i, 0], self.node[i, 1])

    @ti.kernel
    def vector_linear_interpolation(self, f: ti.template(), R: ti.template()):
        """! 定义在二维平面上的向量函数的线性插值
        """
        for i in range(self.node.shape[0]):
            R[i, 0], R[i, 1] = f(self.node[i, 0], self.node[i, 1])

    @ti.kernel
    def surface_linear_scalar_interpolation(self, f: ti.template(), R: ti.template()):
        """! 定义在二维曲面上的标量函数的线性插值
        """
        for i in range(self.node.shape[0]):
            R[i] = f(self.node[i, 0], self.node[i, 1], self.node[i, 2])

    @ti.kernel
    def surface_linear_vector_interpolation(self, f: ti.template(), R: ti.template()):
        """! 定义在二维曲面上的向量函数的线性插值
        """
        for i in range(self.node.shape[0]):
            R[i, 0], R[i, 1] = f(self.node[i, 0], self.node[i, 1], self.node[i, 2])


    @ti.kernel
    def surface_cell_mass_matrix(self, S: ti.template()):
        """
        组装曲面三角形网格上的最低次单元刚度矩阵， 
        这里的曲面是指三角形网格的节点几何维数为 3
        """
        pass
