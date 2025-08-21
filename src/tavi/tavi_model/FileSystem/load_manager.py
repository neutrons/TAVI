import logging
import os
from typing import Dict, Iterable, Optional

import tavi.tavi_model.FileSystem.load_ornl as load_ornl
from tavi.tavi_model.FileSystem.tavi_class import TaviProject

logger = logging.getLogger("TAVI")


class LoadManager:
    """
    TAVI file management system.

    - Normalizes user input (folder and files) and resolves paths.
    - Ranks likely facility using facility-provided scoring functions.
    - Loads data through the chosen facility loader.

    Parameters
    ----------
    data_folder : str | os.PathLike | None
        Folder path containing data files (optional if `data_file` is provided).
    data_file : str | os.PathLike | Iterable[str | os.PathLike] | None
        A single file or an iterable of files.
    facility : str | None
        If provided, bypass ranking and force this facility (must be supported).
    """

    # currently supported facilities
    supported_facilities = ["ORNL", "ANSTO", "ILL", "PSI", "NIST", "MLZ"]

    def __init__(
        self,
        data_folder: Optional[os.PathLike | str] = None,
        data_file: Optional[os.PathLike | str | Iterable[os.PathLike | str]] = None,
        facility: Optional[str] = None,
    ) -> None:
        self.data_folder = data_folder
        self.data_files = data_file
        self.facility = facility

    # TO DO
    def rank_facility(self) -> str:
        """
        Determine and set the most likely facility using plugin score functions.

        Returns
        -------
        str
            The selected facility name.
        """
        scores: Dict[str, float] = {name: 0.0 for name in self.supported_facilities}
        for facility in self.supported_facilities:
            try:
                scores[facility] = load_ornl.score(self.data_folder, self.data_files)
            except Exception as exc:
                logger.exception("Can't compute proper score for %s: %s", facility, exc)
                scores[facility] = float("-inf")

            # TO DO: extend to other facilities
            # elif facility == "ANSTO":
            #     score[facility] = load_ansto.score(self.data_folder, self.datafiles)
            # elif facility == "ILL":
            #     score[facility] = load_ill.score(self.data_folder, self.datafiles)
            # elif facility == "PSI":
            #     score[facility] = load_psi.score(self.data_folder, self.datafiles)
            # elif facility == "NIST":
            #     score[facility] = load_nist.score(self.data_folder, self.datafiles)
            # elif facility == "MLZ":
            #     score[facility] = load_mlz.score(self.data_folder, self.datafiles)
        self.facility = max(scores, key=scores.get)
        return self.facility

    def load(self) -> TaviProject:
        """
        load the either a folder or a single file or a list of files and returns the TaviProject class
        """
        # if self.facility is None, call the ranking system
        if not self.facility:
            self.rank_facility()
        match self.facility:
            case "ORNL":
                return load_ornl.load_ornl(self.data_folder, self.data_files)

            # TO DO: extend to other facilities
            # case "ILL":
            # return load_ill.load(self.data_folder, self.data_files)
            # case "ANSTO":
            # return load_ansto.load(self.data_folder, self.data_files)
            # case "PSI":
            # return load_psi.load(self.data_folder, self.data_files)
            # case "NIST":
            # return load_nist.load(self.data_folder, self.data_files)
            # case "MLZ":
            # return load_ill.load(self.data_folder, self.data_files)
