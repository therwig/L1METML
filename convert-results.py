import uproot
import h5py
import numpy as np

hf = h5py.File("output_MET_DY.h5", "r")
f = uproot.recreate("output_MET_DY.root")

f["Events"] = uproot.newtree( hf )
f["Events"].extend( {k:np.array(hf[k]) for k in hf.keys()} )

f.close()
hf.close()
