
# preprocessimages.m README

## Description: 
* function preprocesses images, then saves the images as .tif and .mat files in 'organizehybs' folder.

## Inputs: 
1. experiment name (used for naming files)
2. experiment directory
	* main directory with HybCycle_[i] images, etc.)
3. position or fov
4. folder array
	* (array for number of folders ex. 0:4 for [0,1,2,3,4] folders.

## Output: 
1. returns processed images as I
	* folder by channel cell array of images
	* saves the images in the project folder as a tiff file and a mat data file
2. Saves dapi aligned images as 'AllHybRegistration.tif'
3. Saves raw images hybIms, dapiIms, dapiTform as imagesHybDapi-pos[0-9]-[experimentName]-[date].mat
4. Saves processed images I, dapiTform, position as preProcessedData-pos[0-9]-[experimentName]-[date].mat

## Dependencies
1. Matlab Version R2018a used
2. Fiji installed
	* can be downloaded at https://imagej.net/Fiji/Downloads

## Requirements
1. images need at least 4 Z-slices for the alignment
2. mij.jar in Fiji.app/jars
	* Download at http://bigwww.epfl.ch/sage/soft/mij/mij.jar
3. 16 GB of memory for number of folders totaling more than 40 
4. Access to save background subtracted images in experimentDir

## Options: 
1. useBackgroundImages: boolean to use background images for subtraction
2. backgroundFolderName: change default folder name in experimentDir from 'final_background' to new folder
3. dapiRefPath: path to align images to certain dapi image. If null, images align to dapi in HybCycle_0 of the experiment.
4. imageJBackSubtract: boolean to use imageJ rolling ball subtraction.
	* automatically set to true
	* useful from leveling background around spots or cells
5. subtractBackground: boolean to subtract background using background images.
	* useful for high background images or autofluorescence.
6. saveProcessedHybIms: boolean to save processed Hyb images
7. divideIms: divides dapi images into 4 to use gradient descent for 3d alignment
	* good for stacks with less than 16 zslices
	* usually don't need to use
0. default Values for optional arguments: [true, 'initial_background', [], false, false, false]; 

## Examples:

#### Preprocess images 
* background images in folder 'initial_background'
* background images used for uneven illumination correction
* imageJBackSubtract is set to true and will level background

```Matlab 
experimentDir = 'I:\2019-06-27_KPC_pool_auto\Experiment';
experimentName = 'KPC_pool-2019-06-27'
position = 0; % positions 0
folderArray = 0:7; % folders 0 - 7 

I = preprocessimages(experimentName, experimentDir, position, folderArray);
```

#### Different Folder Name for Background
* change the folder name of the background images to 'final_background'

```Matlab 
useBackground = true;
backgroundFolder = 'final_background';
I = preprocessimages(experimentName, experimentDir, position, folderArray, useBackground, backgroundFolder);
```

#### No background images
```Matlab 
I = preprocessimages(experimentName, experimentDir, position, folderArray, false);
```

#### Save processed tifs and use Background
* default background folder: 'initial_background' (can set to null)
* And subtract background to get rid of background noise

```Matlab 
experimentDir = 'I:\2019-06-27_KPC_pool_auto\Experiment';
experimentName = 'KPC_pool-2019-06-27';
position = 0;
folderArray = 0:7;
useBackground = true;
subtractBackground = true;
backgroundFolder = []; % default will be used
dapiRefPath = ''; % use HybCycle_0 in experiment directory
imageJBackSubtract = true;
saveProcessIms = true;
I = preprocessimages(experimentName, experimentDir, position, folderArray, ...
	useBackground, backgroundFolder, dapiRefPath, imageJBackSubtract, ...
	subtractBackground, saveProcessIms);
```

#### Change the dapi reference image path
* want to align images to another experiment: experiment2\HybCycle_0
* use initial dapi image in first hyb cycle
* use same settings from previous example with defaults of using the background image for correcting uneven illuminatin bias, imageJ background subtraction, and no background subtraction from background images.

```Matlab
dapiRefPath = 'I:\experiment2\HybCycle_0';
I = preprocessimages(experimentName, experimentDir, position, folderArray, ...
	useBackground, backgroundFolder, dapiRefPath);
```

#### No ImageJ Rolling Ball Background Subtraction
```Matlab
imageJBackSubtract = false;
I = preprocessimages(experimentName, experimentDir, position, folderArray, ...
	useBackground, backgroundFolder, dapiRefPath, imageJBackSubtract);
```

#### Preprocess positions 0,1,2,3: make a for loop
	* Store Images in a 1 by 4 cell array for each position (field of view).
```Matlab 
fovArray = 0:3;
I = cell(1, length(fovArray));
for position = fovArray
	I{position + 1} = preprocessimages(experimentName, experimentDir, position, folderArray, ...
		useBackground); % need to add 1 because Matlab starts indices at 1
end
```
