# dna-seqfish-plus
code generated and used to generate data for seqFISH experiments

## dependencies
1. Fiji
	* download at: https://imagej.net/Fiji/Downloads
	* add the scripts folder to the path
```Matlab
addpath('path\Fiji.app\scripts', '-end');
```
2. bfmatlab (already in io package)
	* download at: https://downloads.openmicroscopy.org/bio-formats/5.3.4/artifacts/bfmatlab.zip
	* add to MATLAB path
	* add "bioformats_package.jar" to javapath
```Matlab
javaaddpath('path\bioformats_package.jar');
```

## packages
| package  | Description |
| ------------- | ------------- |
| align  | dapi or fiducial alignment  |
| io  | TIFF readers and save wrappers  |
| immunofluorescence  | retrieve intensity profile of images  |
| spots  | spot detectors  |
| threshold  | manual and auto  |
| preprocess | preprocessing filters |
| process  | grab spots and organize  |
| decode  | decode genes  |
| segment | io for ilastik masks, assign points to masks  |
| example  | scripts used |
