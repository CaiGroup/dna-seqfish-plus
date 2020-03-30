# process

## seqimprocess

### Description:
* grabs points, segments based on RoiSets.zip, and outputs csv of locations.
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
1. points, pointspercel .mat 
3. .csv output of data for point locations
### Requirements
1. Matlab Version R2019a
2. Processed images
3. Segmentation in 2d (RoiSet.zip) or 3d (labeled images)
4. Threshold 
* matrix of hybcycles by channels

### Optional Arguments
1. sqrtradius:
* square root radius for colocalizing points
* default is 6, which is equal to 2.45
2. typedots:
* type of dot detection: 
	* 'log' is the default LoG filter
	* 'exons' uses laplacian of gaussians to find inflection points
	* 'introns' finds local peaks
3. superresolve:
* type of super resolution
	* 'none' (default) keeps the original points
	* 'gaussian' uses a gausian to find peak
	* 'radial' uses radial center algorithm to find peak 
		* much faster than gaussian
### Examples
#### Output positions for points 
```Matlab 
experimentDir, 'D:\exampleExp\';
experimentName = 'exampleExperiment2';
chArray = 1:2; % use channels 1 and 2
folderArray = 0:19; % number of folders or hybcycles
position = 0;
typedots = 'exons';
superres = 'radial'; 
sqrtradius = 6;
segoption = '2d'; % no 3d option yet
[pointspercell, points] = seqimprocess(experimentDir, experimentName, ...
    processedImages, position, folderArray, chArray, segoption, ...
    sqrtradius, typedots, superres)
```
