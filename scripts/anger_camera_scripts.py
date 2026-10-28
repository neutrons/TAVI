from tavi.library.camera.loader import AngerCameraLoader
from tavi.library.camera.algorigthms import AngerCameraAlgo
import numpy as np

loader = AngerCameraLoader()
loader.parse_folder_xml("/Users/qmc/Neutrons/tasdata/DEMAND/exp940_For_Qianli")
# loader.load_ub("/Users/qmc/Neutrons/tasdata/DEMAND/exp940_For_Qianli/4K_UB_spice.dat")
loader.set_ub("0.02122675 0.00785842 0.07885838 0.19478760 0.00692103 -0.00721790 0.10426915 0.22170969 -0.00256979", ub_name="matlab")
ub = "matlab"
anger_camera_data = loader.data

ac = AngerCameraAlgo(anger_camera_data, wavelength=1.533)

scan_name = "HB3A_exp940_scan0009"
ac.show_scans(scan_name = scan_name, ub_name="matlab", save = "test.pdf",
              peak_kwargs=dict(sigma1 = 8, sigma2 =24, threshold = 1.0, min_counts = 500), 
              vmin = 0, vmax = 1)


# scan_name = "HB3A_exp940_scan0009_0039.xml"
# detector = ac.convert_camera_hkl(scan_name, ub)
# filtered = ac.gaussian_filter_diff(detector=detector, sigma1=6, sigma2=12)
# peaks = ac.find_peaks(detector=detector)
# print(peaks)