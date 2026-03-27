"""Cooper-Nathans Method."""
import numpy as np
from tavi.library.component.collimators import Collimators
from tavi.library.component.monochromater import Monochromator
from tavi.library.component.analyzer import Analyzer

class CooperNathans:
    """
    Cooper-Nathans method. 

    Follows Popvici 1975 paper. Calculate matrix G, F, C, A, B.
    """
    
    NUM_MONO = 1
    NUM_ANA = 1
    IDX_MONO = {"0H": 0, "0V":1}
    IDX_ANA_0_H = {"0H":2, "0V":3}

    # 4 soller slits collimators
    NUM_COLLS = 4
    IDX_COLLS = {"0H":0, "0V":2, "1H":1, "1V":3, "2H":4, "2V":6, "3H":5, "3V":7}
    
    def mat_g(self, collimators : Collimators) -> np.ndarray:
        """Matrix G, 8 x 8 diagonal matrix about the horizontal and vertical divergence of 4 collimators."""
        mat_g = np.zeros(self.NUM_COLLS * 2, self.NUM_COLLS * 2)
        
        mat_g[self.IDX_COLLS["0H"], self.IDX_COLLS["0H"]] = 1.0/collimators.pre_mono_h ** 2
        mat_g[self.IDX_COLLS["0V"], self.IDX_COLLS["0V"]] = 1.0/collimators.pre_mono_v ** 2
        
        mat_g[self.IDX_COLLS["1H"], self.IDX_COLLS["1H"]] = 1.0/collimators.pre_sample_h ** 2
        mat_g[self.IDX_COLLS["1V"], self.IDX_COLLS["1V"]] = 1.0/collimators.pre_sample_v ** 2

        mat_g[self.IDX_COLLS["2H"], self.IDX_COLLS["2H"]] = 1.0/collimators.post_sample_h ** 2
        mat_g[self.IDX_COLLS["2V"], self.IDX_COLLS["2V"]] = 1.0/collimators.post_sample_v ** 2

        mat_g[self.IDX_COLLS["3H"], self.IDX_COLLS["3H"]] = 1.0/collimators.post_ana_h ** 2
        mat_g[self.IDX_COLLS["3V"], self.IDX_COLLS["3V"]] = 1.0/collimators.post_ana_v ** 2

        return mat_g



    def mat_f(self, monochromator : Monochromator, analyzer : Analyzer) -> np.ndarray:
        """
        Matrix F, a 4 x 4 diagonal matrix considering divergence of monochromator and analyzer.
        """

        mat_f = np.zerors((self.NUM_MONO + self.NUM_ANA) * 2, (self.NUM_MONO + self.NUM_ANA) * 2)
        mat_f[self.IDX_MONO["0H"], self.IDX_MONO["0H"]] = 1.0 / monochromator.mosaic_h **2
        mat_f[self.IDX_MONO["0V"], self.IDX_MONO["0V"]] = 1.0 / monochromator.mosaic_v **2
        mat_f[self.IDX_ANA["0H"], self.IDX_ANA["0H"]] = 1.0 / analyzer.mosaic_h **2
        mat_f[self.IDX_ANA["0V"], self.IDX_ANA["0V"]] = 1.0 / analyzer.mosaic_v **2
        return mat_f

    def mat_a(self, ki, kf, theta_m, theta_a):
        """Calculate Matrix A. 6 x 8 matrix, Y = AU. Transform collmiator's divergence to ki-kf."""
        pass
