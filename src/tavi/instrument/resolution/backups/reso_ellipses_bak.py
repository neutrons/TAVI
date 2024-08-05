import numpy as np
import numpy.linalg as la
import mpl_toolkits.mplot3d as mplot3d
import matplotlib
import matplotlib.pyplot as plt
from tavi.utilities import *


np.set_printoptions(floatmode="fixed", precision=4)


class ResoEllipsoid(object):
    """Manage the resolution ellipoid

    Attributs:
        frame (str): "q", "hkl", or "proj"
        STATUS (None | bool): True if resolution calculation is successful
        mat (float): 4 by 4 resolution matrix
        r0 (float | None): Normalization factor

    Methods:

    """

    g_eps = 1e-8

    def __init__(self):

        self.q = None
        self.en = None
        self.frame = None
        self.projection = None

        self.STATUS = None
        self.mat = None
        self.r0 = None
        self.coh_fwhms = None
        self.incoh_fwhms = None
        self.principal_fwhms = None

    # def __init__(self, mat, r0):

    #     self.mat = mat
    #     self.r0 = r0
    #     self.coh_fwhms = None
    #     self.incoh_fwhms = None
    #     self.principal_fwhms = None

    #     if np.isnan(r0) or np.isinf(r0) or np.isnan(mat.any()) or np.isinf(mat.any()):
    #         self.STATUS = False
    #     else:
    #         self.STATUS = True

    def ellipsoid_volume(self):
        """volume of the ellipsoid"""
        pass

    @staticmethod
    def quadric_proj(quadric, idx):
        """projects along one axis of the quadric"""

        if np.abs(quadric[idx, idx]) < ResoEllipsoid.g_eps:
            return np.delete(np.delete(quadric, idx, axis=0), idx, axis=1)

        # row/column along which to perform the orthogonal projection
        vec = 0.5 * (quadric[idx, :] + quadric[:, idx])  # symmetrise if not symmetric
        vec /= np.sqrt(quadric[idx, idx])  # normalise to indexed component
        proj_op = np.outer(vec, vec)  # projection operator
        ortho_proj = quadric - proj_op  # projected quadric

        # comparing with simple projection
        # rank = len(quadric)
        # vec /= np.sqrt(np.dot(vec, vec))
        # proj_op = np.outer(vec, vec)
        # ortho_proj = np.dot((np.identity(rank) - proj_op), quadric)

        # remove row/column that was projected out
        # print("\nProjected row/column %d:\n%s\n->\n%s.\n" % (idx, str(quadric), str(ortho_proj)))
        return np.delete(np.delete(ortho_proj, idx, axis=0), idx, axis=1)

    def coh_fwhms_calc(self):
        """Coherent FWHMs"""
        reso = self.mat
        fwhms = []

        for i in range(len(reso)):
            fwhms.append(sig2fwhm / np.sqrt(reso[i, i]))
        fwhms = np.array(fwhms)
        self.coh_fwhms = fwhms

    def incoh_fwhms_calc(self):
        """Incoherent FWHMs"""

        reso = self.mat
        proj_q_para = ResoEllipsoid.quadric_proj(reso, 3)
        proj_q_para = ResoEllipsoid.quadric_proj(proj_q_para, 2)
        proj_q_para = ResoEllipsoid.quadric_proj(proj_q_para, 1)

        proj_q_perp = ResoEllipsoid.quadric_proj(reso, 3)
        proj_q_perp = ResoEllipsoid.quadric_proj(proj_q_perp, 2)
        proj_q_perp = ResoEllipsoid.quadric_proj(proj_q_perp, 0)

        proj_q_up = ResoEllipsoid.quadric_proj(reso, 3)
        proj_q_up = ResoEllipsoid.quadric_proj(proj_q_up, 1)
        proj_q_up = ResoEllipsoid.quadric_proj(proj_q_up, 0)

        proj_en = ResoEllipsoid.quadric_proj(reso, 2)
        proj_en = ResoEllipsoid.quadric_proj(proj_en, 1)
        proj_en = ResoEllipsoid.quadric_proj(proj_en, 0)

        fwhms = np.array(
            [
                1.0 / np.sqrt(np.abs(proj_q_para[0, 0])) * sig2fwhm,
                1.0 / np.sqrt(np.abs(proj_q_perp[0, 0])) * sig2fwhm,
                1.0 / np.sqrt(np.abs(proj_q_up[0, 0])) * sig2fwhm,
                1.0 / np.sqrt(np.abs(proj_en[0, 0])) * sig2fwhm,
            ]
        )
        self.incoh_fwhms = fwhms

    def principal_fwhms_calc(self):
        """FWHMs of principal axes in 4D"""

        evals, _ = la.eig(self.mat)
        fwhms = np.array(1.0 / np.sqrt(np.abs(evals)) * sig2fwhm)
        self.principal_fwhms = fwhms

    @staticmethod
    def descr_ellipse(quadric):
        """2D ellipse"""

        evals, evecs = la.eig(quadric)  # evecs normalized already
        fwhms = 1.0 / np.sqrt(np.abs(evals)) * sig2fwhm
        # angles = np.array([np.arctan2(evecs[1][0], evecs[0][0])])

        return (fwhms, evecs)

    @staticmethod
    def _gen_ellipse_cut(mat, axes=(0, 1)):
        """generate a 2D ellipsis by removing other axes"""
        q_res = mat
        for i in (3, 2, 1, 0):
            if i not in axes:
                q_res = np.delete(np.delete(q_res, i, axis=0), i, axis=1)

        fwhms, vec = ResoEllipsoid.descr_ellipse(q_res)
        return (fwhms, vec)

    @staticmethod
    def _gen_ellipse_project(mat, axes=(0, 1)):
        """generate a 2D ellipsis by projecting out other axes"""
        q_res = mat
        # Qres_QxE_proj = np.delete(np.delete(self.mat, 2, axis=0), 2, axis=1)
        for i in (3, 2, 1, 0):
            if i not in axes:
                q_res = ResoEllipsoid.quadric_proj(q_res, i)

        fwhms, rot = ResoEllipsoid.descr_ellipse(q_res)
        return (fwhms, rot)

    def calc_ellipses(self, verbose=True):
        """Calculate FWHMs and rotation angles"""
        self.coh_fwhms_calc()
        self.incoh_fwhms_calc()
        self.principal_fwhms_calc()

        if verbose:
            print(f"Coherent-elastic fwhms:{self.coh_fwhms}")
            print(f"Incoherent-elastic fwhms: {self.incoh_fwhms}")
            print(f"Principal axes fwhms: {self.principal_fwhms}")

        # 2d sliced ellipses

        fwhms_QxE, rot_QxE = ResoEllipsoid._gen_ellipse_cut(self.mat, axes=(0, 3))
        fwhms_QyE, rot_QyE = ResoEllipsoid._gen_ellipse_cut(self.mat, axes=(1, 3))
        fwhms_QzE, rot_QzE = ResoEllipsoid._gen_ellipse_cut(self.mat, axes=(2, 3))

        fwhms_QxQy, rot_QxQy = ResoEllipsoid._gen_ellipse_cut(self.mat, axes=(0, 1))
        fwhms_QyQz, rot_QyQz = ResoEllipsoid._gen_ellipse_cut(self.mat, axes=(1, 2))
        fwhms_QxQz, rot_QxQz = ResoEllipsoid._gen_ellipse_cut(self.mat, axes=(0, 2))

        if verbose:
            print(f"2d Qx,E slice fwhms and vectors: {fwhms_QxE}, {rot_QxE}")
            print(f"2d Qy,E slice fwhms and vectors:{fwhms_QyE},{rot_QyE}")
            print(f"2d Qz,E slice fwhms and vectors:{fwhms_QzE},{rot_QzE}")
            print(f"2d Qx,Qy slice fwhms and vectors:{fwhms_QxQy},{rot_QxQy}")
            print(f"2d Qy,Qz slice fwhms and vectors:{fwhms_QyQz},{rot_QyQz}")
            print(f"2d Qx,Qz slice fwhms and vectors:{fwhms_QxQz},{rot_QxQz}")

        # 2d projected ellipses
        # Qres_QxE_proj = np.delete(np.delete(self.mat, 2, axis=0), 2, axis=1)

        (fwhms_QxE_proj, rot_QxE_proj) = ResoEllipsoid._gen_ellipse_project(self.mat, axes=(0, 3))
        (fwhms_QyE_proj, rot_QyE_proj) = ResoEllipsoid._gen_ellipse_project(self.mat, axes=(1, 3))
        (fwhms_QzE_proj, rot_QzE_proj) = ResoEllipsoid._gen_ellipse_project(self.mat, axes=(2, 3))
        (fwhms_QxQy_proj, rot_QxQy_proj) = ResoEllipsoid._gen_ellipse_project(self.mat, axes=(0, 1))
        (fwhms_QyQz_proj, rot_QyQz_proj) = ResoEllipsoid._gen_ellipse_project(self.mat, axes=(1, 2))
        (fwhms_QxQz_proj, rot_QxQz_proj) = ResoEllipsoid._gen_ellipse_project(self.mat, axes=(0, 2))

        if verbose:
            print(f"2d Qx,E projection fwhms and vectors: {fwhms_QxE_proj}, {rot_QxE_proj}")
            print(f"2d Qy,E projection fwhms and vectors: {fwhms_QyE_proj}, {rot_QyE_proj}")
            print(f"2d Qz,E projection fwhms and vectors: {fwhms_QzE_proj}, {rot_QzE_proj}")
            print(f"2d Qx,Qy projection fwhms and vectors: {fwhms_QxQy_proj}, {rot_QxQy_proj}")
            print(f"2d Qy,Qz projection fwhms and vectors: {fwhms_QyQz_proj}, {rot_QyQz_proj}")
            print(f"2d Qx,Qz projection fwhms and vectors: {fwhms_QxQz_proj}, {rot_QxQz_proj}")

        results = {
            "fwhms_QxE": fwhms_QxE,
            "rot_QxE": rot_QxE,
            "fwhms_QyE": fwhms_QyE,
            "rot_QyE": rot_QyE,
            "fwhms_QzE": fwhms_QzE,
            "rot_QzE": rot_QzE,
            "fwhms_QxQy": fwhms_QxQy,
            "rot_QxQy": rot_QxQy,
            "fwhms_QyQz": fwhms_QyQz,
            "rot_QyQz": rot_QyQz,
            "fwhms_QxQz": fwhms_QxQz,
            "rot_QxQz": rot_QxQz,
            "fwhms_QxE_proj": fwhms_QxE_proj,
            "rot_QxE_proj": rot_QxE_proj,
            "fwhms_QyE_proj": fwhms_QyE_proj,
            "rot_QyE_proj": rot_QyE_proj,
            "fwhms_QzE_proj": fwhms_QzE_proj,
            "rot_QzE_proj": rot_QzE_proj,
            "fwhms_QxQy_proj": fwhms_QxQy_proj,
            "rot_QxQy_proj": rot_QxQy_proj,
            "fwhms_QyQz_proj": fwhms_QyQz_proj,
            "rot_QyQz_proj": rot_QyQz_proj,
            "fwhms_QxQz_proj": fwhms_QxQz_proj,
            "rot_QxQz_proj": rot_QxQz_proj,
        }

        return results

    def plot_ellipses(
        self,
        ellis,
        verbose=True,
        plot_results=True,
        file="",
        dpi=600,
        ellipse_points=128,
        use_tex=False,
    ):
        """Plot 2D and 3D ellipses"""
        matplotlib.rc("text", usetex=use_tex)

        def ellfkt(rad, vecs, phi):
            return np.dot(
                vecs,
                np.array(
                    [
                        rad[0] * np.cos(phi),
                        rad[1] * np.sin(phi),
                    ],
                ),
            )

        phi = np.linspace(0, 2.0 * np.pi, ellipse_points)

        ell_QxE = ellfkt(ellis["fwhms_QxE"] * 0.5, ellis["rot_QxE"], phi)
        ell_QyE = ellfkt(ellis["fwhms_QyE"] * 0.5, ellis["rot_QyE"], phi)
        ell_QzE = ellfkt(ellis["fwhms_QzE"] * 0.5, ellis["rot_QzE"], phi)
        ell_QxQy = ellfkt(ellis["fwhms_QxQy"] * 0.5, ellis["rot_QxQy"], phi)
        ell_QyQz = ellfkt(ellis["fwhms_QyQz"] * 0.5, ellis["rot_QyQz"], phi)
        ell_QxQz = ellfkt(ellis["fwhms_QxQz"] * 0.5, ellis["rot_QxQz"], phi)

        ell_QxE_proj = ellfkt(ellis["fwhms_QxE_proj"] * 0.5, ellis["rot_QxE_proj"], phi)
        ell_QyE_proj = ellfkt(ellis["fwhms_QyE_proj"] * 0.5, ellis["rot_QyE_proj"], phi)
        ell_QzE_proj = ellfkt(ellis["fwhms_QzE_proj"] * 0.5, ellis["rot_QzE_proj"], phi)
        ell_QxQy_proj = ellfkt(ellis["fwhms_QxQy_proj"] * 0.5, ellis["rot_QxQy_proj"], phi)
        ell_QyQz_proj = ellfkt(ellis["fwhms_QyQz_proj"] * 0.5, ellis["rot_QyQz_proj"], phi)
        ell_QxQz_proj = ellfkt(ellis["fwhms_QxQz_proj"] * 0.5, ellis["rot_QxQz_proj"], phi)

        labelQpara = "Qpara (1/A)"
        labelQperp = "Qperp (1/A)"
        labelQup = "Qup (1/A)"

        if use_tex:
            labelQpara = "$Q_{\parallel}$ (\AA$^{-1}$)"
            labelQperp = "$Q_{\perp}$ (\AA$^{-1}$)"
            labelQup = "$Q_{up}$ (\AA$^{-1}$)"

        # Qpara, E axis

        fig = plt.figure()
        subplot_QxE = fig.add_subplot(231)
        subplot_QxE.set_xlabel(labelQpara)
        subplot_QxE.set_ylabel("E (meV)")
        subplot_QxE.plot(ell_QxE[0], ell_QxE[1], c="black", linestyle="dashed")
        subplot_QxE.plot(ell_QxE_proj[0], ell_QxE_proj[1], c="black", linestyle="solid")
        subplot_QxE.grid(alpha=0.6)

        # Qperp, E axis
        subplot_QyE = fig.add_subplot(232)
        subplot_QyE.set_xlabel(labelQperp)
        subplot_QyE.set_ylabel("E (meV)")
        subplot_QyE.plot(ell_QyE[0], ell_QyE[1], c="black", linestyle="dashed")
        subplot_QyE.plot(ell_QyE_proj[0], ell_QyE_proj[1], c="black", linestyle="solid")
        subplot_QyE.grid(alpha=0.6)
        # Qup, E axis
        subplot_QzE = fig.add_subplot(233)
        subplot_QzE.set_xlabel(labelQup)
        subplot_QzE.set_ylabel("E (meV)")
        subplot_QzE.plot(ell_QzE[0], ell_QzE[1], c="black", linestyle="dashed")
        subplot_QzE.plot(ell_QzE_proj[0], ell_QzE_proj[1], c="black", linestyle="solid")
        subplot_QzE.grid(alpha=0.6)
        # Qpara, Qperp axis
        subplot_QxQy = fig.add_subplot(234)
        subplot_QxQy.set_xlabel(labelQpara)
        subplot_QxQy.set_ylabel(labelQperp)
        subplot_QxQy.plot(ell_QxQy[0], ell_QxQy[1], c="black", linestyle="dashed")
        subplot_QxQy.plot(ell_QxQy_proj[0], ell_QxQy_proj[1], c="black", linestyle="solid")
        subplot_QxQy.grid(alpha=0.6)
        subplot_QxQy.set_aspect("equal")
        # Qperp, Qup axis
        subplot_QyQz = fig.add_subplot(235)
        subplot_QyQz.set_xlabel(labelQperp)
        subplot_QyQz.set_ylabel(labelQup)
        subplot_QyQz.plot(ell_QyQz[0], ell_QyQz[1], c="black", linestyle="dashed")
        subplot_QyQz.plot(ell_QyQz_proj[0], ell_QyQz_proj[1], c="black", linestyle="solid")
        subplot_QyQz.grid(alpha=0.6)
        subplot_QyQz.set_aspect("equal")
        # Qpara, Qup axis
        subplot_QxQz = fig.add_subplot(236)
        subplot_QxQz.set_xlabel(labelQpara)
        subplot_QxQz.set_ylabel(labelQup)
        subplot_QxQz.plot(ell_QxQz[0], ell_QxQz[1], c="black", linestyle="dashed")
        subplot_QxQz.plot(ell_QxQz_proj[0], ell_QxQz_proj[1], c="black", linestyle="solid")
        subplot_QxQz.grid(alpha=0.6)
        subplot_QxQz.set_aspect("equal")

        plt.tight_layout()

        # 3d plot
        fig3d = plt.figure()
        subplot3d = fig3d.add_subplot(111, projection="3d")

        subplot3d.set_xlabel(labelQpara)
        subplot3d.set_ylabel(labelQperp)
        # subplot3d.set_ylabel(labelQup)
        subplot3d.set_zlabel("E (meV)")

        # xE
        subplot3d.plot(ell_QxE[0], ell_QxE[1], zs=0.0, zdir="y", c="black", linestyle="dashed")
        subplot3d.plot(ell_QxE_proj[0], ell_QxE_proj[1], zs=0.0, zdir="y", c="black", linestyle="solid")
        # yE
        subplot3d.plot(ell_QyE[0], ell_QyE[1], zs=0.0, zdir="x", c="black", linestyle="dashed")
        subplot3d.plot(ell_QyE_proj[0], ell_QyE_proj[1], zs=0.0, zdir="x", c="black", linestyle="solid")
        # zE
        # subplot3d.plot(ell_QzE[0], ell_QzE[1], zs=0., zdir="x", c="black", linestyle="dashed")
        # subplot3d.plot(ell_QzE_proj[0], ell_QzE_proj[1], zs=0., zdir="x", c="black", linestyle="solid")
        # xy
        subplot3d.plot(ell_QxQy[0], ell_QxQy[1], zs=0.0, zdir="z", c="black", linestyle="dashed")
        subplot3d.plot(ell_QxQy_proj[0], ell_QxQy_proj[1], zs=0.0, zdir="z", c="black", linestyle="solid")

        # save plots to files
        if file != "":
            import os

            splitext = os.path.splitext(file)
            file3d = splitext[0] + "_3d" + splitext[1]

            if verbose:
                print('Saving 2d plot to "%s".' % file)
                print('Saving 3d plot to "%s".' % file3d)
            fig.savefig(file, dpi=dpi)
            fig3d.savefig(file3d, dpi=dpi)

        # show plots
        if plot_results:
            plt.show()
