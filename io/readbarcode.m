function barcodekey = readbarcode(barcodePath, varargin)
% barcodekey reads an excel file and coverts the data into a structure of
% barcodekey with values of barcode matrix and cell array of gene names.
%
% Inputs: (barcode excel file path, optional parameters 'header' and 'no
% header' to indicate if header is present. Default is 'no header'.
%
% Outputs: barcodekey structure
%
% Author: Nico Pierson
% Date: November 12, 2018

%% Set up optional Parameters for z-slice index
    numvarargs = length(varargin);
    if numvarargs > 1
        error('streamline:readbarcode:TooManyInputs', ...
            'requires at most 1 optional inputs');
    end
    
    % Error for type of arguments
    if numvarargs >= 1 
        if ~ischar(varargin{1}) 
            error('streamline:readbarcode:WrongInput', ...
                'paramter input requires type string');
        elseif ~strcmp(varargin{1}, 'header') && ~strcmp(varargin{1}, 'no header') 
            error('streamline:readbarcode:WrongInput', ...
                'parameter input requires type string: "header" or "no header"');
        end
    end

    % set defaults for optional inputs
    optargs = {'no header'};
    
    % now put these defaults into the valuesToUse cell array, 
    % and overwrite the ones specified in varargin.
    optargs(1:numvarargs) = varargin;
    
    % Default Value of ref image is 1
    [headerPresent] = optargs{:};

    %% Import Data from Excel File and Store in Structure
    f = importdata(barcodePath);
    barcodekey.barcode = f.data;
    numberOfGenes = size(f.textdata, 1);

    switch headerPresent % Switch statement to see if header is present
        case 'header'
            barcodekey.names = f.textdata(2:numberOfGenes, 1);
        case 'no header'
            barcodekey.names = f.textdata(1:numberOfGenes, 1);
        otherwise
            error(['Invalid optional argument, ', varargin{1}]);
    end
    
end