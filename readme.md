# dna-seqfish-plus
This repostory contains the scripts used in processing the images and barcode calling for the dna seqFISH experiment and the source data which consists of the processed data. 
March 27, 2020

## Getting Started
* Download all the contents of the folder and add it to your Matlab path.

### Prerequisites
* MATLAB 2019

### Dependencies
1. radialcenter.m by Dr. Parthasarathy, which can be downloaded [here](https://media.nature.com/original/nature-assets/nmeth/journal/v9/n7/extref/nmeth.2071-S2.zip).
2. Fiji
	* download at:
	* add the scripts folder to the path
```Matlab
addpath('path\Fiji.app\scripts', '-end');
```
3. bfmatlab
	* download at: https://downloads.openmicroscopy.org/bio-formats/5.3.4/artifacts/bfmatlab.zip
	* add to MATLAB path
	* add "bioformats_package.jar" to javapath
```Matlab
javaaddpath('path\bioformats_package.jar');
```

## Running the Code
*Overview of the steps for the DNA seqFISH pipeline*
### Prerequisites
* Region of interests (ROI) holds polygon cell segmentations created from ImageJ software
* thresholding values after Preprocessing Step

## packages
example scripts can be found in each package
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
| extra  | extra scripts used for sequential thresholding and spot finding |

## License
Free for non-commercial and academic research. The software is bound by the licensing rules of California Institute of Technology (Caltech)

## Image Processing Team
* Yodai Takei
* Sheel Shah - decoding, spot finding, and 3D radial center
* Nico Pierson - writing, cleaning, debugging, and validating code
* Jonathan White - fiduciary alignment


## Contact
* Contact the corresponding author: Long Cai (lcai@caltech.edu) for any inquiry.



