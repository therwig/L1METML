import uproot
import h5py
import numpy as np

hf = h5py.File("output_MET_DY.h5", "r")
f = uproot.recreate("output_MET_DY.root")

f["Events"] = uproot.newtree()

for k in hf.keys():
    print 'saving a branch',k
    f["Events"].extend({k:np.array(hf[k])})

f.close()
hf.close()
