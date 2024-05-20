    @staticmethod
    def convert_spice_to_hdf5(path_to_spice_folder, path_to_hdf5):
        """Load data from spice folder.

        Args:
            path_to_spice_folder (str): spice folder, ends with '/'
            path_to_hdf5 (str): path to hdf5 data file, ends with '.h5'
        """
        exp_info = [
            "experiment",
            "experiment_number",
            "proposal",
            "users",
            "local_contact",
        ]

        p = Path(path_to_spice_folder)

        with h5py.File(path_to_hdf5, "w") as f:

            scans = sorted((p / "Datafiles").glob("*"))
            instrument_str, exp_str = scans[0].parts[-1].split("_")[0:2]

            # read in exp_info from the first scan and save as attibutes of the file
            _, _, headers, _ = TAVI_Data._read_spice(scans[0])
            ipts = headers["proposal"]
            exp_id = "IPTS" + ipts + "_" + instrument_str

            grp_data = f.create_group("data_" + exp_id)
            grp_processed_data = f.create_group("processed_data")
            grp_fit = f.create_group("fits")
            grp_plot = f.create_group("plots")

            for k, v in headers.items():
                if k in exp_info:
                    grp_data.attrs[k] = v

            # read scans into dataset1
            for scan in scans:  # ignoring unused keys
                spice_data, col_headers, headers, unused = TAVI_Data._read_spice(scan)

                scan_num = ((scan.parts[-1].split("_"))[-1]).split(".")[0]
                scan_id = exp_str + "_" + scan_num
                scan_entry = grp_data.create_group(scan_id)
                scan_entry.attrs["scan_id"] = scan_id

                for k, v in headers.items():
                    if k not in exp_info:  # ignore common keys in single scans
                        if "," in v and k != "scan_title":  # vectors
                            scan_entry.attrs[k] = np.array([float(v0) for v0 in v.split(",")])
                        elif v.replace(".", "").isnumeric():  # numebrs only
                            if v.isdigit():  # int
                                scan_entry.attrs[k] = int(v)
                            else:  # float
                                scan_entry.attrs[k] = float(v)
                        # separate COM/FWHM and its errorbar
                        elif k == "Center of Mass":
                            com, e_com = v.split("+/-")
                            scan_entry.attrs["COM"] = float(com)
                            scan_entry.attrs["COM_err"] = float(e_com)
                        elif k == "Full Width Half-Maximum":
                            fwhm, e_fwhm = v.split("+/-")
                            scan_entry.attrs["FWHM"] = float(fwhm)
                            scan_entry.attrs["FWHM_err"] = float(e_fwhm)
                        else:  # other crap, keep as is
                            if k not in exp_info:
                                scan_entry.attrs[k] = v

                if spice_data.ndim == 1:  # empty data or 1 point only
                    if len(spice_data):  # 1 point only
                        for idx, col_header in enumerate(col_headers):
                            scan_entry.create_dataset(col_header, data=spice_data[idx])
                    else:  # empty
                        pass
                else:  # nomarl data
                    for idx, col_header in enumerate(col_headers):
                        scan_entry.create_dataset(col_header, data=spice_data[:, idx])
