# threshold
Package used to manually threshold images
Date: 8/28/2019

## thresholdbych

### Description: 
* grabs images from the experimentDir, or uses input images to manually threshold; outputs threshold as matrices in .mat files

### Inputs: 
1. experiment name (used for naming files)
2. experiment directory
	* main directory with HybCycle_[i] images, etc.)
3. position
4. folder array
	* (array for number of folders ex. 0:4 for [0,1,2,3,4] folders.
5. typedots:
	* type of dot detection: 
		* 'log' (default) log filter in 3d that uses a laplacian of gaussians
		* 'exons' uses laplacian of gaussians to find inflectionpoints
		* 'introns' finds local peaks
### Optional Arguments
1. images:
	* use processed or raw images in workspace
### Output: 
1. threshold for each channel
	* saved in 'threshold' folder in experiment dir

### Dependencies
1. Matlab Version R2018a
2. Fiji installed
	* can be downloaded at https://imagej.net/Fiji/Downloads

### Requirements: Preprocessing
1. images need at least 4 Z-slices for the alignment
2. mij.jar in Fiji.app/jars
	* Download at http://bigwww.epfl.ch/sage/soft/mij/mij.jar
3. 16 GB of memory maximum 40 folders 
4. Access to save background subtracted images in experimentDir
### Examples
#### Loading images from experiment
```Matlab 
experimentDir = 'D:\exampleUser\exampleDir\experiment1';
experimentName = 'experiment1';
position = 1;
folderArray = 0:19;
typedots = 'exons';
threshold = thresholdbych(experimentDir, experimentName, position, ...
    folderArray, typedots)
```
#### Using images from workspace
```Matlab 
experimentDir = 'D:\exampleUser\exampleDir\experiment1';
experimentName = 'experiment1';
position = 1;
folderArray = 0:19;
typedots = 'exons'; 
threshold = thresholdbych(experimentDir, experimentName, position, ...
    folderArray, typedots, images)
```
