# dna-seqfish-plus
* code generated and used to generate data for seqFISH experiments

## dependencies
1. Fiji
	* download at:https://imagej.net/Fiji/Downloads
	* add the scripts folder to the path
```Matlab
addpath('path\Fiji.app\scripts', '-end');
```
2. bfmatlab
	* download at: https://downloads.openmicroscopy.org/bio-formats/5.3.4/artifacts/bfmatlab.zip
	* add to MATLAB path
	* add "bioformats_package.jar" to javapath
```Matlab
javaaddpath('path\bioformats_package.jar');
```

## packages
1. align
* functions to align images for dapi and fiducial markers
2. process
* functions to process images
3. decoding
* decode points
4. segment
* segment cells
5. immunofluorescence
* retrieve pixel intensities of antibody stains
6. spots
* dot detection
7. threshold
* manual and auto threshold
8. preprocess
* preprocess images: imagej background subtraction, uneven illumination correction
9. examples
* example scripts
10. io
* input/output TIFF, csv, and mat files
