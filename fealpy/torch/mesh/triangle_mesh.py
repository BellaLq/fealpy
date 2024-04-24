
from typing import Optional

import numpy as np
import torch
from torch import Tensor

from .. import logger
from . import functional as F
from . import mesh_kernel as K
from .mesh_base import HomoMeshDataStructure, HomoMesh

_dtype = torch.dtype
_device = torch.device


class TriangleMeshDataStructure(HomoMeshDataStructure):
    def __init__(self, NN: int, cell: Tensor):
        super().__init__(NN, 2, cell)
        # constant tensors
        kwargs = {'dtype': cell.dtype, 'device': cell.device}
        self.localEdge = torch.tensor([(1, 2), (2, 0), (0, 1)], **kwargs)
        self.localFace = torch.tensor([(1, 2), (2, 0), (0, 1)], **kwargs)
        self.ccw = torch.tensor([0, 1, 2], **kwargs)

        self.localCell = torch.tensor([
            (0, 1, 2),
            (1, 2, 0),
            (2, 0, 1)], **kwargs)

        self.construct()

    def total_face(self):
        return self.cell[..., self.localFace].reshape(-1, 2)

    def construct(self) -> None:
        kwargs = {'dtype': self.itype, 'device': self.device}
        NC = self.cell.shape[0]
        NFC = self.cell.shape[1]

        totalFace = self.total_face()
        _, i0_np, j_np = np.unique(
            torch.sort(totalFace, dim=1)[0].cpu().numpy(),
            return_index=True,
            return_inverse=True,
            axis=0
        )
        self.face = totalFace[i0_np, :] # this also adds the edge in 2-d meshes
        NF = i0_np.shape[0]

        i1_np = np.zeros(NF, **kwargs)
        i1_np[j_np] = np.arange(NFC*NC, **kwargs)

        self.cell2edge = torch.from_numpy(j_np).to(self.device).reshape(NC, NFC)
        self.cell2face = self.cell2edge

        face2cell_np = np.stack([i0_np//NFC, i1_np//NFC, i0_np%NFC, i1_np%NFC], axis=-1)
        self.face2cell = torch.from_numpy(face2cell_np).to(self.device)
        self.edge2cell = self.face2cell

        logger.info(f"Mesh toplogy relation constructed, with {NF} edge (or face), "
                    f"on device {self.device}")


class TriangleMesh(HomoMesh):

    def shape_function(self, bc: Tensor, p: int=1, mi: Optional[Tensor]=None):
        TD = bc.shape[-1] - 1
        mi = mi or F.multi_index_matrix(p, TD, dtype=self.ds.itype, device=self.device)
        return K.simplex_shape_function(bc, p, mi)