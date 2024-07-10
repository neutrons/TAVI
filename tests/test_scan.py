import numpy as np
from tavi.tavi_data.scans import Scan


def test_plot_gen():
    s2 = Scan()
    s2.scan_info = dict(def_x="s1", def_y="detector")
    s2.data = dict(s1=np.arange(0, 1, 0.2), detector=np.arange(100, 300, 40))
    x, y, xerr, yerr, _, _, _ = s2.generate_curve()
    print(x)
    assert np.allclose(x, np.arange(0, 1, 0.2))
    # assert y == np.arange(100, 300, 40)
    # assert xerr == 0
    # assert yerr == 0


if __name__ == "__main__":
    test_plot_gen()
