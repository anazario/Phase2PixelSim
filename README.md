# Phase2PixelSim
Software for Phase 2 Pixel Simulation 

## Set Up Code
```
cmsrel CMSSW_9_3_2
cd CMSSW_9_3_2/src/
cmsenv
git clone git@github.com:anazario/Phase2PixelSim.git
scram b -j9
```

## Get Geometry XML Files
```
mkdir VFPix
cd VFPix
git init 
git remote add -f origin git@github.com:OSU-CMS/VFPix.git
git config core.sparsecheckout true
echo MonteCarlo/data/OT613_200_IT4025_opt6s3l/ >> .git/info/sparse-checkout
echo MonteCarlo/data/OT613_200_IT4025_opt7s4l/ >> .git/info/sparse-checkout
echo MonteCarlo/data/OT613_200_IT4025_opt8s3l/ >> .git/info/sparse-checkout
git pull origin master
```

## Test interactively with each geometry

..*Program: Phase2PixelSim/Phase2PixelStubs/python/Phase2PixelStubs_cfg.py
..*using options: Phase2PixelStubs_cfg.py start=(int) finish=(int) geometry=(string)

### Available Options 
start: number of file to start on (default 0)
finish: number of file to end on (default 1)
geometry: geometry to run over (default opt8s4l, can also run opt6s3l, opt7s4l and opt8s3l)
Additional geometries can be added in the geomDict.py. They are added as a list of step 3 files
to be run over using the batchList.py script.
### Example
```
cd Phase2PixelSim/Phase2PixelStubs/python
cmsRun Phase2PixelStubs_cfg.py start=0 finish=10 geometry=opt7s4l
```

## Run on Condor
Geometries are run separately
Program: Phase2PixelSim/Phase2PixelStubs/python/condor/condorSubmit.py
options: condorSubmit.py -n (int), -g (string), -c (no argument)
-n: number of files per job (default 5)
-g: geometry to be run by the jobs (only ones recognized: opt8s4l (default), opt6s3l, opt7s4l and opt8s3l)
-c: create the condor_submit.txt 

### Run the code
```
cd Phase2PixelSim/Phase2PixelStubs/python/condor
python condorSubmit.py (default options)
```
