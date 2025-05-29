import matplotlib.pyplot as plt

from tavi.data.scan import Scan

scan006 = Scan.from_spice("./scripts/test_scan/", scan_num=6)
scan006.plot()

plt.show()
