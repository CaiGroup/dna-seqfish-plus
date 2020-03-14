# immuno

## avgintpercell

### Description: 
* function avgintpercell.m takes the segmentations, bins the segmented cells into 2x2 grid and calculates the average intensity across the positions, folder, channels, and cells.

### Inputs: 
1. experiment name (used for naming files)
2. experiment directory
	* main directory with HybCycle_[i] images, etc.)
3. Raw or processed images
	* folder by channel cell array
4. position
5. folder array
	* (array for number of folders ex. 0:4 for [0,1,2,3,4] folders.
6. channel array
	* 1:2 for channels 1 and 2
7. segmentaiton option:
	* '2d' or '3d'

### Output: 
1. .csv file of average intensity in 'analysis' folder
2. .mat file of segmentation and boundaries for each position

### Dependencies
1. Matlab Version R2018a
2. Fiji installed
	* can be downloaded at https://imagej.net/Fiji/Downloads
3. Segmentation in 2d (RoiSet.zip) or 3d (labeled images)

### Requirements: Preprocessing
1. images need at least 4 Z-slices for the alignment
2. mij.jar in Fiji.app/jars
	* Download at http://bigwww.epfl.ch/sage/soft/mij/mij.jar
3. 16 GB of memory maximum 40 folders 
4. Access to save background subtracted images in experimentDir
 

### Examples:

#### First, Preprocess images for Position 0
* don't use background images nor imageJ rolling ball background subtraction
```Matlab
experimentDir = 'I:\2019-07-25-E14-DNA-seqFISH+rep2-2-DNAFISH-plate2';
experimentName = 'ImmunoFluorescence-2019-07-25-E14-DNA-seqFISH+rep2-2-DNAFISH-plate2';
position = 0;
folderArray = 0:20;
imageJBackSubtract = false;
useBackground = false;
backgroundFolder = [];
dapiRefPath = 'I:\2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped\HybCycle_0'; % want another dapi as the reference image to align
[I, ~, ~] = preprocessimages(experimentName, experimentDir, ...
position, folderArray, useBackground, backgroundFolder, dapiRefPath, imageJBackSubtract);
```

#### Output Average Intensity .csv file 
* use same variables from preprocessing
```Matlab 
chArray = 1:2;
folderArray = 0:19;
segoption = '2d';
avgintpercell(experimentDir, experimentName, I, position, folderArray, chArray, segoption); % ADD OPTION IF IMAGE IS SHIFTED
```