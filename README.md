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

* Configuration file location: ```Phase2PixelSim/Phase2PixelStubs/python/Phase2PixelStubs_cfg.py```

### Available Options 
* start: number of file to start on (default 0).
* finish: number of file to end on (default 1).
* geometry: geometry to run over (default opt8s4l). 
### Using options 

```cmsRun Phase2PixelStubs_cfg.py start=(int) finish=(int) geometry=(string)```

### Available Geometries
* opt8s4l: Standard Geometry 4025, 8 small discs and 4 large discs.
* opt8s3l: 8 small discs and 3 large discs.
* opt7s4l: 7 small discs and 4 large discs.
* opt6s3l: 6 small discs and 3 large discs.

### Example
```
cd Phase2PixelSim/Phase2PixelStubs/python
cmsRun Phase2PixelStubs_cfg.py start=0 finish=10 geometry=opt7s4l
```

## Run on Condor
**Geometries are run one at a time**
* Program: ```Phase2PixelSim/Phase2PixelStubs/python/condor/condorSubmit.py```
* options: condorSubmit.py -n (int), -g (string), -c (no argument).

  * -n: number of files per job (default 5).
  * -g: geometry to be run by the jobs (default opt8s4l). See available ones above.
  * -c: create the condor_submit.txt submit file without running. 

### Condor Example
```
cmsenv
cd $CMSSW_BASE/Phase2PixelSim/Phase2PixelStubs/python/condor
python condorSubmit.py -n 10 -g opt8s3l 
```

## Create Plots
The file for producing the plots is ```Phase2PixelSim/Phase2PixelStubs/python/Plot3.C``` and uses the class definitions in ```Plot2.h```. Currently the ```Plot3.C``` script makes very specific plots comparing the standard ('opt8s4l') geometry with all of the other geometries ('opt8s3l', 'opt7s4l' and 'opt6s3l'). Plots include total stub distribution, number of stubs per eta, and number of stubs in each of the barrel layers, endcap double discs and in each of the endcap disc rings for both the positive and negative regions of the CMS detector.

1. hadd the condor files for each seperate geometry.
2. Modify the lines in ```Plot3.C``` corresponding to the location of the hadded files and each of the ROOT file names.
3. Run the script:
 ```root -l Plot3.C```

## Adding new geometries
Geometries from MC simulation can be added in the form of a txt file containing the appropriate "step 3" root files.  

Additional geometries can be added in: 
```Phase2PixelSim/Phase2PixelStubs/python/geomDict.py```

### Create text file from finished condor jobs:
 ```
 cmsenv
 cd $CMSSW_BASE/Phase2PixelSim/Phase2PixelStubs/python
 python batchList.py -d <condor_jobs_pathname> -l
 mv <text_filename>.txt GeomRootFiles/
 ```
 
 The location of the text file can be added to the dictionary in geomDict.py as:
 
 ```dict['<geometry option name>'] = environ["CMSSW_BASE"]+"/src/Phase2PixelSim/Phase2PixelStubs/python/GeomRootFiles/<filename>.txt"```
 
 ### Example
 ```dict['opt6s3l'] = environ["CMSSW_BASE"]+"/src/Phase2PixelSim/Phase2PixelStubs/python/GeomRootFiles/OT613_200_IT4025_opt6s3l_step3.txt"```
