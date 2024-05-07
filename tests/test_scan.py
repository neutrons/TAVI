import numpy as np
from tavi.scans import Scan


def test_get_metadata():
    s1 = Scan()
    s1.metadata = dict(IPTS=1234, exp=567, def_x="s1")
    # print(s1)
    # assert s1.metadata["def_x"] == "test"
    md = s1.get_metadata()
    assert md["def_x"] == "s1"


def test_plot_gen():
    s2 = Scan()
    s2.metadata = dict(def_x="s1", def_y="detector")
    s2.data = dict(s1=np.arange(0, 1, 0.2), detector=np.arange(100, 300, 40))
    x, y, xerr, yerr = s2.plot_gen()
    print(x)
    assert np.allclose(x, np.arange(0, 1, 0.2))
    # assert y == np.arange(100, 300, 40)
    # assert xerr == 0
    # assert yerr == 0
