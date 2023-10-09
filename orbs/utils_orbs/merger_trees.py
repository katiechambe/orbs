import numpy as np
from utils_orbs.orb_paths import SetupPaths


class TraceMergerTree:

    def __init__(
        self, 
        snapshot,
        subfindID,
        sim = "TNG",
        physics ="hydro",
        **kwargs
        ):
        """
        Identifies and pulls merger tree for a single subhalo

        Parameters
        ----------
        snapshot: int
            the number of the snapshot with the corresponding subhalo ID
        subfindID: int
            the ID number of the subhalo at the corresponding snapshot
        sim: str
            "Illustris" or "TNG"
            to specify which simulation
        physics: str
            "dark" or "hydro"
            to specify which simulation
        kwargs: dict
            little_h: h varies for each simulation!
        """

        SetupPaths.__init__(self)

        self.snapshot = snapshot
        self.subfindID = subfindID
        self.sim = sim
        self.physics = physics
        self.kwargs = kwargs
        self.little_h = self.kwargs.pop("little_h", 0.704)

        # defining the simulation path from paths.py
        if self.sim == "Illustris":
            totalSnaps = 135
            from utils.readtreeHDF5Py3 import TreeDB
            if self.physics == "dark":
                self.treepath = self.path_illustrisdark_trees
            elif self.physics == "hydro":
                self.treepath = self.path_illustrishydro_trees
                
        elif self.sim == "TNG":
            totalSnaps = 99
            from utils.readtreeHDF5_public import TreeDB
            if self.physics == "dark":
                self.treepath = self.path_tngdark_trees
            elif self.physics == "hydro":
                self.treepath = self.path_tnghydro_trees

        treeDirectory = self.treepath

        tree = TreeDB(treeDirectory)
        pastbranch = tree.get_main_branch( 
            self.snapshot, 
            self.subfindID
            # keysel=['SnapNum', 'SubhaloMass', 'SubhaloPos', 'SubhaloVel', 'SubhaloID', 'SubfindID']
            )
        futurebranch = tree.get_future_branch(
                               self.snapshot,
                               self.subfindID)

        self.pastbranch = pastbranch
        self.futurebranch = futurebranch
        
        self.pastkeys = np.array(list(self.pastbranch.__dict__.keys()))
        self.futurekeys = np.array(list(self.futurebranch.__dict__.keys()))
        
        self.mergedbranch = {}
        for key in self.pastkeys[np.isin(self.pastkeys,self.futurekeys)]:
        #print(type(tree1.futurebranch.__getattribute__(key)))
                self.mergedbranch[key] = np.concatenate([self.futurebranch.__getattribute__(key)[:-1],self.pastbranch.__getattribute__(key)])

        
#         self.snaphist = np.concatenate([futurebranch.SnapNum[0:-1],pastbranch.SnapNum])
#         self.missingsnaps = np.arange(0,totalSnaps+1,1)[~np.isin(np.arange(0,totalSnaps+1,1),self.snaphist)]
        
        
        
#         self.snaps = branch.SnapNum
#         self.masses = branch.SubhaloMass
#         self.positions = branch.SubhaloPos # note this is in comoving!
#         self.velocities = branch.SubhaloVel
#         self.id = branch.SubhaloID 
#         self.subfindIDTree = branch.SubfindID

#         self.masses_phys = self.masses / self.little_h
        
    @property
    def maxmass(self):
        """
        Max mass of the subhalo
        -- note: this only considers current and previous snapshots --

        Parameters:
        -----------
        None

        Outputs:
        --------
        maxmass: float
            the maximum mass previously achieved by a subhalo
        maxsnap: int
            the snapshot at which max mass occurs 
        maxredshift: float
            the maximum redshift at which max mass occurs
        """
        maxmass = max(self.masses_phys)
        maxmass_mask =  max(self.masses_phys)==self.masses_phys
        maxsnap = self.snaps[maxmass_mask][0]
        return maxmass, maxsnap


