
from ..nodetype import CNodeType, PortConf, DataType
from .utils import get_mesh_class

__all__ = [
    "Box2d",
    "Box3d"
]


class Box2d(CNodeType):
    r"""Create a mesh in a box-shaped 2D area.

    Inputs:
        mesh_type (str): Type of mesh to granerate.
        domain (tuple[float, float, float, float], optional): Domain.
        nx (int, optional): Segments on x direction.
        ny (int, optional): Segments on y direction.

    Outputs:
        mesh (MeshType): The mesh object created.
    """
    TITLE: str = "矩形网格"
    PATH: str = "网格.构造"
    INPUT_SLOTS = [
        PortConf("mesh_type", DataType.MENU, 0, title="网格类型", default="triangle", items=["triangle", "quadrangle"]),
        PortConf("domain", DataType.NONE, title="区域"),
        PortConf("nx", DataType.INT, title="X 分段数", default=10, min_val=1),
        PortConf("ny", DataType.INT, title="Y 分段数", default=10, min_val=1)
    ]
    OUTPUT_SLOTS = [
        PortConf("mesh", DataType.MESH, title="网格")
    ]

    @staticmethod
    def run(mesh_type, domain, nx, ny):
        MeshClass = get_mesh_class(mesh_type)
        kwds = {"nx": nx, "ny": ny}
        if domain is not None:
            kwds["box"] = domain
        return MeshClass.from_box(**kwds)


class Box3d(CNodeType):
    r"""Create a mesh in a box-shaped 3D area.

    Inputs:
        mesh_type (str): Type of mesh to granerate.
        domain (tuple[float, float, float, float, float, float], optional): Domain.
        nx (int, optional): Segments on x direction.
        ny (int, optional): Segments on y direction.
        nz (int, optional): Segments on z direction.

    Outputs:
        mesh (MeshType): The mesh object created.
    """
    TITLE: str = "长方体网格"
    PATH: str = "网格.构造"
    INPUT_SLOTS = [
        PortConf("mesh_type", DataType.MENU, 0, title="网格类型", default="tetrahedron", items=["tetrahedron", "hexahedron"]),
        PortConf("domain", DataType.NONE, title="区域"),
        PortConf("nx", DataType.INT, title="X 分段数", default=10, min_val=1),
        PortConf("ny", DataType.INT, title="Y 分段数", default=10, min_val=1),
        PortConf("nz", DataType.INT, title="Z 分段数", default=10, min_val=1)
    ]
    OUTPUT_SLOTS = [
        PortConf("mesh", DataType.MESH, title="网格")
    ]

    @staticmethod
    def run(mesh_type, domain, nx, ny, nz):
        MeshClass = get_mesh_class(mesh_type)
        kwds = {"nx": nx, "ny": ny, "nz": nz}
        if domain is not None:
            kwds["box"] = domain
        return MeshClass.from_box(**kwds)
