# L1METML

Setup:
```
cmsrel CMSSW_10_6_1_patch2
cd CMSSW_10_6_1_patch2/src
cmsenv
git clone https://github.com/jmduarte/L1METML
cd L1METML
```


Pre-process root file input to calculate quantities wrt Z boson for training:
```
g++ reprocess_inputs.cc `root-config --cflags` `root-config --glibs`
./a.out
```

Convert root file to h5:
```
python convert-uproot.py
```

Train with h5 inputs:
```
python train.py
```

Convert back to root for evaluation:
```
python convert-results.py
```

Run pyRoot analysis scripts:
```
python quickana.py
```

