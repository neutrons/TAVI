import numpy as np
import numpy.linalg as la
import mpl_toolkits.mplot3d as mplot3d
import matplotlib
import matplotlib.pyplot as plt
from tavi.utilities import *


np.set_printoptions(floatmode="fixed", precision=4)


class Reso_Ellipsoid(object):
    """Manage the resolution ellipoid

    Attributs:
        STATUS (bool): True if resolution calculation is successful
        mat (float): 4 by 4 resolution matrix
        r0 (float | None): Normalization factor

    Methods:

    """

    g_eps = 1e-8

    def __init__(self, mat, r0):

        self.mat = mat
        self.r0 = r0
        self.coh_fwhms = None
        self.incoh_fwhms = None
        self.principal_fwhms = None

        if np.isnan(r0) or np.isinf(r0) or np.isnan(mat.any()) or np.isinf(mat.any()):
            self.STATUS = False
        else:
            self.STATUS = True

    def ellipsoid_volume(self):
        """volume of the ellipsoid"""
        pass

    @staticmethod
    def quadric_proj(quadric, idx):
        """projects along one axis of the quadric"""

        if np.abs(quadric[idx, idx]) < Reso_Ellipsoid.g_eps:
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
        proj_q_para = Reso_Ellipsoid.quadric_proj(reso, 3)
        proj_q_para = Reso_Ellipsoid.quadric_proj(proj_q_para, 2)
        proj_q_para = Reso_Ellipsoid.quadric_proj(proj_q_para, 1)

        proj_q_perp = Reso_Ellipsoid.quadric_proj(reso, 3)
        proj_q_perp = Reso_Ellipsoid.quadric_proj(proj_q_perp, 2)
        proj_q_perp = Reso_Ellipsoid.quadric_proj(proj_q_perp, 0)

        proj_q_up = Reso_Ellipsoid.quadric_proj(reso, 3)
        proj_q_up = Reso_Ellipsoid.quadric_proj(proj_q_up, 1)
        proj_q_up = Reso_Ellipsoid.quadric_proj(proj_q_up, 0)

        proj_en = Reso_Ellipsoid.quadric_proj(reso, 2)
        proj_en = Reso_Ellipsoid.quadric_proj(proj_en, 1)
        proj_en = Reso_Ellipsoid.quadric_proj(proj_en, 0)

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

        evals, evecs = la.eig(quadric)
        fwhms = 1.0 / np.sqrt(np.abs(evals)) * sig2fwhm
        angles = np.array([np.arctan2(evecs[1][0], evecs[0][0])])

        return (fwhms, angles * rad2deg, evecs)

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
        Qres_QxE = np.delete(np.delete(self.mat, 2, axis=0), 2, axis=1)
        Qres_QxE = np.delete(np.delete(Qres_QxE, 1, axis=0), 1, axis=1)
        [fwhms_QxE, angles_QxE, rot_QxE] = Reso_Ellipsoid.descr_ellipse(Qres_QxE)

        Qres_QyE = np.delete(np.delete(self.mat, 2, axis=0), 2, axis=1)
        Qres_QyE = np.delete(np.delete(Qres_QyE, 0, axis=0), 0, axis=1)
        [fwhms_QyE, angles_QyE, rot_QyE] = Reso_Ellipsoid.descr_ellipse(Qres_QyE)

        Qres_QzE = np.delete(np.delete(self.mat, 1, axis=0), 1, axis=1)
        Qres_QzE = np.delete(np.delete(Qres_QzE, 0, axis=0), 0, axis=1)
        [fwhms_QzE, angles_QzE, rot_QzE] = Reso_Ellipsoid.descr_ellipse(Qres_QzE)

        Qres_QxQy = np.delete(np.delete(self.mat, 3, axis=0), 3, axis=1)
        Qres_QxQy = np.delete(np.delete(Qres_QxQy, 2, axis=0), 2, axis=1)
        [fwhms_QxQy, angles_QxQy, rot_QxQy] = Reso_Ellipsoid.descr_ellipse(Qres_QxQy)

        if verbose:
            print(f"2d Qx,E slice fwhms and slope angle: {fwhms_QxE}, {angles_QxE[0]}")
            print(f"2d Qy,E slice fwhms and slope angle:{fwhms_QyE},{angles_QyE[0]}")
            print(f"2d Qz,E slice fwhms and slope angle:{fwhms_QzE},{angles_QzE[0]}")
            print(f"2d Qx,Qy slice fwhms and slope angle:{fwhms_QxQy},{angles_QxQy[0]}")

        # 2d projected ellipses
        Qres_QxE_proj = np.delete(np.delete(self.mat, 2, axis=0), 2, axis=1)
        Qres_QxE_proj = Reso_Ellipsoid.quadric_proj(Qres_QxE_proj, 1)
        [fwhms_QxE_proj, angles_QxE_proj, rot_QxE_proj] = Reso_Ellipsoid.descr_ellipse(Qres_QxE_proj)

        Qres_QyE_proj = np.delete(np.delete(self.mat, 2, axis=0), 2, axis=1)
        Qres_QyE_proj = Reso_Ellipsoid.quadric_proj(Qres_QyE_proj, 0)
        [fwhms_QyE_proj, angles_QyE_proj, rot_QyE_proj] = Reso_Ellipsoid.descr_ellipse(Qres_QyE_proj)

        Qres_QzE_proj = np.delete(np.delete(self.mat, 1, axis=0), 1, axis=1)
        Qres_QzE_proj = Reso_Ellipsoid.quadric_proj(Qres_QzE_proj, 0)
        [fwhms_QzE_proj, angles_QzE_proj, rot_QzE_proj] = Reso_Ellipsoid.descr_ellipse(Qres_QzE_proj)

        Qres_QxQy_proj = Reso_Ellipsoid.quadric_proj(self.mat, 3)
        Qres_QxQy_proj = np.delete(np.delete(Qres_QxQy_proj, 2, axis=0), 2, axis=1)
        [fwhms_QxQy_proj, angles_QxQy_proj, rot_QxQy_proj] = Reso_Ellipsoid.descr_ellipse(Qres_QxQy_proj)

        if verbose:
            print(f"2d Qx,E projection fwhms and slope angle: {fwhms_QxE_proj}, {angles_QxE_proj[0]}")
            print(f"2d Qy,E projection fwhms and slope angle: {fwhms_QyE_proj}, {angles_QyE_proj[0]}")
            print(f"2d Qz,E projection fwhms and slope angle: {fwhms_QzE_proj}, {angles_QzE_proj[0]}")
            print(f"2d Qx,Qy projection fwhms and slope angle: {fwhms_QxQy_proj}, {angles_QxQy_proj[0]}")

        results = {
            "fwhms_QxE": fwhms_QxE,
            "rot_QxE": rot_QxE,
            "fwhms_QyE": fwhms_QyE,
            "rot_QyE": rot_QyE,
            "fwhms_QzE": fwhms_QzE,
            "rot_QzE": rot_QzE,
            "fwhms_QxQy": fwhms_QxQy,
            "rot_QxQy": rot_QxQy,
            "fwhms_QxE_proj": fwhms_QxE_proj,
            "rot_QxE_proj": rot_QxE_proj,
            "fwhms_QyE_proj": fwhms_QyE_proj,
            "rot_QyE_proj": rot_QyE_proj,
            "fwhms_QzE_proj": fwhms_QzE_proj,
            "rot_QzE_proj": rot_QzE_proj,
            "fwhms_QxQy_proj": fwhms_QxQy_proj,
            "rot_QxQy_proj": rot_QxQy_proj,
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

        ellfkt = lambda rad, rot, phi: np.dot(rot, np.array([rad[0] * np.cos(phi), rad[1] * np.sin(phi)]))

        phi = np.linspace(0, 2.0 * np.pi, ellipse_points)

        ell_QxE = ellfkt(ellis["fwhms_QxE"] * 0.5, ellis["rot_QxE"], phi)
        ell_QyE = ellfkt(ellis["fwhms_QyE"] * 0.5, ellis["rot_QyE"], phi)
        ell_QzE = ellfkt(ellis["fwhms_QzE"] * 0.5, ellis["rot_QzE"], phi)
        ell_QxQy = ellfkt(ellis["fwhms_QxQy"] * 0.5, ellis["rot_QxQy"], phi)

        ell_QxE_proj = ellfkt(ellis["fwhms_QxE_proj"] * 0.5, ellis["rot_QxE_proj"], phi)
        ell_QyE_proj = ellfkt(ellis["fwhms_QyE_proj"] * 0.5, ellis["rot_QyE_proj"], phi)
        ell_QzE_proj = ellfkt(ellis["fwhms_QzE_proj"] * 0.5, ellis["rot_QzE_proj"], phi)
        ell_QxQy_proj = ellfkt(ellis["fwhms_QxQy_proj"] * 0.5, ellis["rot_QxQy_proj"], phi)

        labelQpara = "Qpara (1/A)"
        labelQperp = "Qperp (1/A)"
        labelQup = "Qup (1/A)"

        if use_tex:
            labelQpara = "$Q_{\parallel}$ (\AA$^{-1}$)"
            labelQperp = "$Q_{\perp}$ (\AA$^{-1}$)"
            labelQup = "$Q_{up}$ (\AA$^{-1}$)"

        # Qpara, E axis

        fig = plt.figure()
        subplot_QxE = fig.add_subplot(221)
        subplot_QxE.set_xlabel(labelQpara)
        subplot_QxE.set_ylabel("E (meV)")
        subplot_QxE.plot(ell_QxE[0], ell_QxE[1], c="black", linestyle="dashed")
        subplot_QxE.plot(ell_QxE_proj[0], ell_QxE_proj[1], c="black", linestyle="solid")

        # Qperp, E axis
        subplot_QyE = fig.add_subplot(222)
        subplot_QyE.set_xlabel(labelQperp)
        subplot_QyE.set_ylabel("E (meV)")
        subplot_QyE.plot(ell_QyE[0], ell_QyE[1], c="black", linestyle="dashed")
        subplot_QyE.plot(ell_QyE_proj[0], ell_QyE_proj[1], c="black", linestyle="solid")

        # Qup, E axis
        subplot_QzE = fig.add_subplot(223)
        subplot_QzE.set_xlabel(labelQup)
        subplot_QzE.set_ylabel("E (meV)")
        subplot_QzE.plot(ell_QzE[0], ell_QzE[1], c="black", linestyle="dashed")
        subplot_QzE.plot(ell_QzE_proj[0], ell_QzE_proj[1], c="black", linestyle="solid")

        # Qpara, Qperp axis
        subplot_QxQy = fig.add_subplot(224)
        subplot_QxQy.set_xlabel(labelQpara)
        subplot_QxQy.set_ylabel(labelQperp)
        subplot_QxQy.plot(ell_QxQy[0], ell_QxQy[1], c="black", linestyle="dashed")
        subplot_QxQy.plot(ell_QxQy_proj[0], ell_QxQy_proj[1], c="black", linestyle="solid")
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
