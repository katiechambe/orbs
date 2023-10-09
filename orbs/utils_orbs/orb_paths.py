""" 
Paths for orbits project

Dependencies:
-------------
None
"""

__author__ = "Katie Chamberlain"
__date__   = "October 2023"

class SetupPaths:
    """ 
    Defines the paths that are read into each file
    """

    def __init__(
        self, 
        basedir="/xdisk/gbesla/katiechambe/"
        ):

        self.path_basedir = basedir
        self.path_home = self.path_basedir # + "katie/"
        self.path_orbs = self.path_home + "orbs/"

        # directories for illustris data
        self.path_illustris = self.path_basedir + "Illustris/"
        self.path_illustristng = self.path_basedir + "IllustrisTNG/"
        # group catalog paths
        self.path_illustrisdark = self.path_illustris + "GroupCatalogsDark/"
        self.path_illustrishydro = self.path_illustris + "GroupCatalogsHydro/"
        self.path_tngdark = self.path_illustristng + "TNG100-1-Dark/"
        self.path_tnghydro = self.path_illustristng + "TNG100-1/"
        # merger trees
        self.path_illustrisdark_trees = self.path_illustris + "Illustris-1-Dark-MergerTree/"
        self.path_illustrishydro_trees = self.path_illustris + "Illustris-1-MergerTree/"
        self.path_tngdark_trees = self.path_tngdark + "postprocessing/"
        self.path_tnghydro_trees = self.path_tnghydro + "postprocessing/"
        # matched catalogs
        self.path_tngmatch_V = self.path_illustristng + "TNG100-Matched-V/"
        self.path_tngmatch_N = self.path_illustristng + "TNG100-Matched-Nelson/subhalo_matching_to_dark.hdf5"
        # merger history files
        self.path_tnghydro_mergerhist = self.path_tnghydro + "mergerhistory/"
      

        # orbs directories
        self.path_data = self.path_orbs + "data/"
        self.path_plots = self.path_orbs + "plots/"
        self.path_paper = self.path_orbs + "paper/"