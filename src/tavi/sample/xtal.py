# from typing import Optional

# import numpy as np

# from tavi.sample import Sample
# from tavi.utilities import Peak, UBConf


# class Xtal(Sample):
#     """
#     Singel crystal class

#     Attibutes:
#         type (str): "crystal"
#         ub_peaks : peaks used to determine UB matrix
#         ub_matrix (np.adarray): UB matrix
#         inv_ub_matrix (np.ndarray): inverse of UB matrix
#         in_plane_ref: in plane vector in Qsample frame, goniometers at zero
#         plane_normal: normal vector in Qsample frame, goniometers at zero
#         u (tuple): u vector, (h, k, l) along the beam direction when all goniometer angles are zero
#         v (tuple): v vector, (h ,k, l) in the scattering plane
#     Methods:


#     """

#     def __init__(
#         self,
#         lattice_params=(1, 1, 1, 90, 90, 90),
#     ) -> None:
#         super().__init__(lattice_params)
#         self.type = "crystal"

#         self.ub_peaks: Optional[tuple[Peak, ...]] = None
#         self.u_mat: Optional[np.ndarray] = None
#         self.ub_mat: Optional[np.ndarray] = None
#         # self.inv_ub_mat: Optional[np.ndarray] = None

#         self.plane_normal: Optional[np.ndarray] = None
#         self.in_plane_ref: Optional[np.ndarray] = None

#     @classmethod
#     def from_json(cls, path_to_json):
#         xtal = super().from_json(path_to_json)

#         ub_matrix = xtal.json_dict.get("ub_matrix")
#         if ub_matrix is not None:
#             xtal.ub_mat = np.array(ub_matrix).reshape(3, 3)

#         plane_normal = xtal.json_dict.get("plane_normal")
#         if plane_normal is not None:
#             xtal.plane_normal = np.array(plane_normal)

#         in_plane_ref = xtal.json_dict.get("in_plane_ref")
#         if in_plane_ref is not None:
#             xtal.in_plane_ref = np.array(in_plane_ref)

#         return xtal

#     # @property
#     # def u(self) -> np.ndarray:
#     #     """u vector, in reciprocal lattice unit, along beam"""
#     #     return self.ub_matrix_to_uv(self.ub_mat)[0]

#     # @property
#     # def v(self) -> np.ndarray:
#     #     """
#     #     v vector, in reciprocal lattice unit,
#     #     in the horizaontal scattering plane
#     #     """
#     #     return self.ub_matrix_to_uv(self.ub_mat)[1]

#     def set_orientation(self, ubconf: UBConf) -> None:
#         "Set crystal orientation from UB conf"

#         for key, val in ubconf._asdict().items():
#             if val is not None:
#                 setattr(self, key, val)
