"""Cooper-Nathans Method."""
from sympy import Symbol, Matrix, tan, cos, sin, pprint, pi, det,sqrt
import numpy as np
from tavi.library.experiment.utilities import rotation_matrix_2d

class CooperNathans:
    """
    Cooper-Nathans method.

    Follows Popvici 1975 paper. Calculate matrix G, F, C, A, B.
    """

    NUM_MONO = 1
    NUM_ANA = 1
    IDX_MONO = {"0H": 0, "0V": 1}
    IDX_ANA = {"0H": 2, "0V": 3}

    # 4 soller slits collimators
    NUM_COLLS = 4
    IDX_COLLS = {"0H": 0, "0V": 2, "1H": 1, "1V": 3, "2H": 4, "2V": 6, "3H": 5, "3V": 7}


    def mat_g(self) -> np.ndarray:
        """Matrix G, 8 x 8 diagonal matrix about the horizontal and vertical divergence of 4 collimators."""
        
        mat_g = Matrix.zeros(self.NUM_COLLS * 2, self.NUM_COLLS * 2)
        pre_mono_h = Symbol("pre-mono-h")
        pre_mono_v = Symbol("pre-mono-v")
        mat_g[self.IDX_COLLS["0H"], self.IDX_COLLS["0H"]] = 1.0 / pre_mono_h**2
        mat_g[self.IDX_COLLS["0V"], self.IDX_COLLS["0V"]] = 1.0 / pre_mono_v**2

        pre_sample_h = Symbol("pre-sample-h")
        pre_sample_v = Symbol("pre-sample-v")
        mat_g[self.IDX_COLLS["1H"], self.IDX_COLLS["1H"]] = 1.0 / pre_sample_h**2
        mat_g[self.IDX_COLLS["1V"], self.IDX_COLLS["1V"]] = 1.0 / pre_sample_v**2

        post_sample_h = Symbol("post-sample-h")
        post_sample_v = Symbol("post-sample-v")
        mat_g[self.IDX_COLLS["2H"], self.IDX_COLLS["2H"]] = 1.0 / post_sample_h**2
        mat_g[self.IDX_COLLS["2V"], self.IDX_COLLS["2V"]] = 1.0 / post_sample_v**2

        post_ana_h = Symbol("post-ana-h")
        post_ana_v = Symbol("post-ana-v")
        mat_g[self.IDX_COLLS["3H"], self.IDX_COLLS["3H"]] = 1.0 / post_ana_h**2
        mat_g[self.IDX_COLLS["3V"], self.IDX_COLLS["3V"]] = 1.0 / post_ana_v**2

        return mat_g

    def mat_f(self) -> np.ndarray:
        """
        Matrix F.

        A 4 x 4 diagonal matrix considering divergence of monochromator and analyzer.
        """
        mat_f = Matrix.zeros((self.NUM_MONO + self.NUM_ANA) * 2, (self.NUM_MONO + self.NUM_ANA) * 2)
        mono_mosaic_h = Symbol("mono_mosaic_h")
        mono_mosaic_v = Symbol("mono_mosaic_v")
        mat_f[self.IDX_MONO["0H"], self.IDX_MONO["0H"]] = 1.0 / mono_mosaic_h**2
        mat_f[self.IDX_MONO["0V"], self.IDX_MONO["0V"]] = 1.0 / mono_mosaic_v**2
        
        ana_mosaic_h = Symbol("ana_mosaic_h")
        ana_mosaic_v = Symbol("ana_mosaic_v")
        mat_f[self.IDX_ANA["0H"], self.IDX_ANA["0H"]] = 1.0 / ana_mosaic_h**2
        mat_f[self.IDX_ANA["0V"], self.IDX_ANA["0V"]] = 1.0 / ana_mosaic_v**2
        return mat_f

    def mat_a(self) -> np.ndarray:
        """Calculate Matrix A. 6 x 8 matrix, Y = AU. Transform collmiator's divergence to ki-kf."""
        theta_m = Symbol("theta_m")
        theta_a = Symbol("theta_a")
        ki = Symbol("ki")
        kf = Symbol("kf")
        mat_a = Matrix.zeros(6, 2 * self.NUM_COLLS)
        mat_a[0, self.IDX_COLLS["0H"]] = 0.5 * ki / tan(theta_m)
        mat_a[0, self.IDX_COLLS["1H"]] = -0.5 * ki / tan(theta_m)
        mat_a[1, self.IDX_COLLS["1H"]] = ki
        mat_a[2, self.IDX_COLLS["1V"]] = ki

        mat_a[3, self.IDX_COLLS["2H"]] = 0.5 * kf / tan(theta_a)
        mat_a[3, self.IDX_COLLS["3H"]] = -0.5 * kf / tan(theta_a)
        mat_a[4, self.IDX_COLLS["2H"]] = kf
        mat_a[5, self.IDX_COLLS["2V"]] = kf
        return mat_a

    def mat_b(self) -> np.ndarray:
        """Calculate Matrix B. 4 x 6 matrix, X = BY. Transform from ki-kf to q-frame."""
        phi = Symbol("phi")
        two_theta = Symbol("2theta")
        ki = Symbol("ki")
        kf = Symbol("kf")
        hbar = Symbol("hbar")
        m_n = Symbol("m_n")
        mat_b = Matrix.zeros(4, 6)
        mat_b[0:3, 0:3] = Matrix([[cos(phi), sin(phi), 0], [-sin(phi), cos(phi), 0], [0, 0, 1]])
        mat_b[0:3, 3:6] = Matrix([[cos(phi-two_theta), sin(phi-two_theta), 0], [-sin(phi-two_theta), cos(phi-two_theta), 0], [0, 0, 1]]) * (-1)
        mat_b[3, 0] = ki * hbar / m_n
        mat_b[3, 3] = -kf * hbar / m_n
        return mat_b

    def mat_c(self) -> np.ndarray:
        """Matrix C. 4 x 8 matrix. Constrinat between mono/ana mosaic and collimator divergence."""
        theta_m = Symbol("theta_m")
        theta_a = Symbol("theta_a")
        mat_c = Matrix.zeros((self.NUM_MONO + self.NUM_ANA) * 2, self.NUM_COLLS * 2)
        mat_c[self.IDX_MONO["0H"], self.IDX_COLLS["0H"]] = 0.5
        mat_c[self.IDX_MONO["0H"], self.IDX_COLLS["1H"]] = 0.5
        mat_c[self.IDX_MONO["0V"], self.IDX_COLLS["0V"]] = 0.5 / sin(theta_m)
        mat_c[self.IDX_MONO["0V"], self.IDX_COLLS["1V"]] = -0.5 / sin(theta_m)

        mat_c[self.IDX_ANA["0H"], self.IDX_COLLS["2H"]] = 0.5
        mat_c[self.IDX_ANA["0H"], self.IDX_COLLS["3H"]] = 0.5
        mat_c[self.IDX_ANA["0V"], self.IDX_COLLS["2V"]] = 0.5 / sin(theta_a)
        mat_c[self.IDX_ANA["0V"], self.IDX_COLLS["3V"]] = -0.5 / sin(theta_a)
        return mat_c

    def calculate_at_hkle(self):
        """Calculate the resolution matrix and R0 factor in the local Q frame"""

        mat_h = self.mat_g() + self.mat_c().T @ self.mat_f() @ self.mat_c() 
        mat_h_inv = mat_h.inv()
        mat_ba = self.mat_b() @ self.mat_a()
        # pprint(mat_ba)
        mat_cov = mat_ba @ mat_h_inv @ mat_ba.T
        return mat_cov

if __name__ == "__main__":
    cp = CooperNathans()
    # print("det(G)=")
    # detG = det(cp.mat_g())

    # print("det(F)=")  
    # detF = det(cp.mat_f())

    # print("det(H)=")
    # detH = det(cp.mat_g() + cp.mat_c().T @ cp.mat_f() @ cp.mat_c())

    # pprint((detG*detF))
    # pprint(detH/detG/detF)
    # mat_h = cp.mat_g() + cp.mat_c().T @ cp.mat_f() @ cp.mat_c()
    # pprint(mat_h @ cp.mat_a().inv())
    # print(mat_h.shape)
    # print(mat_ab.shape)
    # pprint(mat_ab.T @ mat_h @ mat_ab)
    # print("det(C.T @ F @ C=")
    # pprint(det(cp.mat_c()@cp.mat_c().T))

    print("C=")
    pprint(cp.mat_c())
    # print("A=")
    # pprint(cp.mat_a())

    # print("B=")
    # pprint(cp.mat_b())

    # mat_cov = cp.calculate_at_hkle()
    # pprint(mat_cov.inv())
    # 4 x 4
    #[x x 0 x]
    #[x x 0 x]
    #[0 0 x 0]
    #[x x 0 x]
    # pprint(mat_cov[0,3])

    # import scipy.constants as co
    # ksq2E = (co.hbar / co.e) ** 2.0 * co.e / 2.0 / co.neutron_mass * 1e23    
    # print(co.hbar**2/co.m_n*1e23/co.e)

