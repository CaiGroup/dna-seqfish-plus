function processImage = imagejbackgroundsubtraction(image, varargin)
% use the imageJ Rolling Ball Background Subtraction on the image, using a
% default of 3 pixels.
%
% To Do List: 
% 1. Add options for stacks of images.
% 2. Delete temporary files for background subtraction...
% - tried fclose, clearvars for image, and MIJ.exit: nothing has worked,
% tried temporary path
%
% Author: Nico Pierson
% Email: nicogpt@caltech.edu
% Date: 1/29/2019

    %% Set up optional Parameters
    
    numvarargs = length(varargin);
    if numvarargs > 4
        error('myfuns:getallbeads:TooManyInputs', ...
            'requires at most 4 optional inputs');
    end

    % set defaults for optional inputs
    optargs = {'tempProcessed', [], 3, false}; % default of using 7 x 7 pixel grid for gaussian function
    
    % now put these defaults into the valuesToUse cell array, 
    % and overwrite the ones specified in varargin.
    optargs(1:numvarargs) = varargin;
    
    % Place optional args in memorable variable names
    [uniqueString, fijiDirectory, pixelRadius, sliding] = optargs{:};

    %% Declare Variables
    % instead of using fijiDirectory, use temporary directory
    if isempty(fijiDirectory)
        fijiDirectory = pwd;
    end
    %tempDirectory = tempdir; % temp direcotyr
    
    % Start Fiji in the background
    Miji(false); 
    
    %% Create image in ImageJ
    namesh = ['C' num2str(1) '-'  num2str(1) '.tif'];
    MIJ.createImage(namesh, image, true);
    % add .tif to uniqueString
    uniqueString = [uniqueString '.tif'];
    saveBackSubImPath = fullfile(fijiDirectory, uniqueString);

    %% Subtract Background using Rolling Ball Algorithm in ImageJ
    try
        if sliding
            MIJ.run('Subtract Background...', ['rolling=' num2str(pixelRadius) ' sliding stack']); % Rolling Background Subtraction ImageJ
        else
            MIJ.run('Subtract Background...', ['rolling=' num2str(pixelRadius) ' stack']); % Rolling Background Subtraction ImageJ
        end
        
        MIJ.run('Save', ['save=[' saveBackSubImPath ']']); % use a unique string to temporarily save image
        MIJ.run('Close All')
    catch
        error('MIJ exited incorrectly: most likely caused by out of memory in the java heap\n');
    end
    
    tempPosition = 0;
    [processImageTemp, ~, ~, ~, ~] = grabimseries(fullfile(fijiDirectory, uniqueString), tempPosition);
    processImage = processImageTemp{1};
    % try to delete temp file
    warning('off','all');
    delete(fullfile(fijiDirectory, uniqueString));
    warning('on','all');

%{
    To do: make function support multiple channels
    
    error in bfmatlab on Alignel's desktop; use grabimseries instead
    
    Warning: Invalid file or directory
'C:\Users\alignell\Desktop\Mike\Fiji.app\scripts\..\plugins\loci_tools.jar'. 
> In javaclasspath>local_validate_dynamic_path (line 271)
  In javaclasspath>local_javapath (line 187)
  In javaclasspath (line 124)
  In javaaddpath (line 71)
  In bfopen (line 103)
  In grabim (line 106)
  In imagejbackgroundsubtraction (line 55)
  In preprocessimages (line 184) 
Warning: Objects of ij/ImagePlus class exist - not clearing java 
> In javaclasspath>local_javapath (line 195)
  In javaclasspath (line 124)
  In javaaddpath (line 71)
  In bfopen (line 103)
  In grabim (line 106)
  In imagejbackgroundsubtraction (line 55)
  In preprocessimages (line 184) 
    
%}