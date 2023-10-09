import numpy as np
import h5py
import sys
import os

class FindHistory:
    """
    Note: this script has only been tested with TNG100-1 run
    ---
    inputs:
    """
    
    
    def __init__(self, historydir):
        self.historydir = historydir
        
#     def __del__(self):
 
        
    def get_info(self, snapnum, subfind_id):
        # check for existence of merger history info:
        if not os.path.exists(f"{self.historydir}MergerHistory_{str(snapnum).zfill(3)}.hdf5"):
            print(f"Path not found:{self.historydir}MergerHistory_{str(snapnum).zfill(3)}.hdf5")
#             sys.exit()
            
        f = h5py.File(f"{self.historydir}MergerHistory_{str(snapnum).zfill(3)}.hdf5", 'r')
        
        datablock = {}
        for key, val in f.items():
            if key == "Header":
                continue
            datablock[key] = np.array(val)[subfind_id]
        
        f.close()
        
        return datablock

#         # Check that a few files/paths exist
#         for rel_path in ['%s.0.hdf5' % (name), 'offsets']:
#             if not os.path.exists(treedir + '/' + rel_path):
#                 print('Path not found: ' + treedir + '/' + rel_path)
#                 sys.exit()

# ,mhfile = h5py.File(f"{paths.path_tnghydro_mergerhist}MergerHistory_040.hdf5","r")