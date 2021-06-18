function outputFile = octExplorerXMLToMatlab(dataRootDir, varargin)
% Converts Explorer XML output to Matlab
%
%  Synopsis
%     outputFile = octExplorerXMLToMatlab(dataRootDir, varargin)
%
% Inputs
%    dataRootDir - Directory containing the XML file
%
% Optional key/val pairs
%    infileSuffix
%    outfileSuffix
%    subjectsToProcess
%
% Outputs
%    outputFile - Filename written to disk
%
% Description
%  Geoff Aguirre collected OCT data using a Heidelberg Spectralis system. A
%  horizontal and vertical macular scan was collected for the left and right
%  eye for each subject. The data were exported from the OCT system and
%  saved in '.E2E' format (which is proprietary to Heidelberg) and in '.vol'
%  format.
% 
%  The '.vol' files were copied to a separate location and then processed by
%  an operator (Kara Cloud) using OCT Explorer v5.0 (on a Macintosh). This
%  analysis yields (for each '.vol' file) a file named
%  '_Surfaces_Retina-JEI-Final.xml'.
% 
%  At Stanford, we also ran OCT Explorer v5.0 on a Mac and wrote out the
%  files.  But the user interactions was surely different.
%
%  The routine 'xmlConversionWrapper.m' converts the '.xml' file to a '.mat'
%  file. The wrapper makes use of the routine 'xml2volmask.m' which is
%  present within 'octSupport', and was written by Jin Gahm from the LONI
%  group, USC.
% 
%  The resulting '.mat' file has the dimensions 1536x97x496, corresponding
%  to the vertical, axial (depth) and horizontal diimensions. Each voxel is
%  given an integer value from 0 - 11, corresponding to a retinal layer
%  (with a value of zero indicating that the voxel does not reside within
%  the retina). The depth dimension is in units of mm, while the transverse
%  (horizontal and vertical) dimensions are in units of degrees of visual
%  angle; our macular acquisitions were 30° wide.
%
% See also
%    ophvlfeat (repository), ophthalmology (directory in scitranApps)
%

%{
% Assuming retinaTOMEAnalysis and ophthalmology are on your path
% Calls this function to convert the Explorer XML outputs to mat-file
  datadir = fullfile(ophRootPath,'iowa','zeissimg');
  octExplorerXMLToMatlab(datadir)
%}
%{
foo = load('P73304206_Macular Cube 512x128_8-19-2020_13-28-54_OS_sn211046_cube_raw_Surfaces_Retina-JEI-Final');
sz = size(foo.mask);
mrvNewGraphWin;
for ii=1:3:sz(2)
  imagesc(double(squeeze(foo.mask(:,ii,:)))); drawnow; pause(0.1); 
end
%}
%{
for ii=1:2:sz(1)
  imagesc(double(squeeze(foo.mask(ii,:,:)))); drawnow; pause(0.1); 
end
%}
%{
for ii=1:10:sz(3)
  imagesc(double(squeeze(foo.mask(:,:,ii)))); drawnow; pause(0.1); 
end
%}
%{
% Show an isosurface surface
% Below I set the 0 values to NaN.  Just testing now.
  mask = foo.mask;
  sz = size(mask);
  mask = double(mask);
  [X,Y,Z] = meshgrid(1:sz(2),1:sz(1),1:sz(3));
  mrvNewGraphWin;
  p = patch(isosurface(X,Y,Z,mask,1.2));
  N = isonormals(X,Y,Z,mask,p);

  % Trying to get into FW.  see t_meshFibersOBJ.m for a successful example
  % Not sure why this doesn't work on the Flywheel side yet.
  FV.vertices = p.Vertices;
  FV.faces    = p.Faces;
  OBJ = objFVN(FV,N);
  disp(FV)
  disp(OBJ)
  fname = fullfile(vistaRootPath,'local','OCT.obj'); 
  objWrite(OBJ,fname);

  p.FaceColor = 'yellow'; p.EdgeColor = 'none';
  daspect([1 1 1]); view(3); axis tight; camlight; lighting gouraud
%}
%{
  foo.mask(foo.mask==0) = NaN;
%}


%% Parse vargin for options passed here

p = inputParser;

% Required
p.addRequired('dataRootDir',@ischar);

% Optional analysis params
p.addParameter('inFileSuffix','_Surfaces_Retina-JEI-Final.xml',@ischar);
p.addParameter('outFileSuffix','_Surfaces_Retina-JEI-Final.mat',@ischar);
p.addParameter('subjectsToProcess',{},@(x)(isempty(x) || iscell(x)));


p.parse(dataRootDir, varargin{:});

%% Convert OCT seg XML to .mat volumes


% Obtain the paths to all of the pupil data files within the specified
% directory, including within sub-drectories.
fileListStruct=subdir(fullfile(dataRootDir,['*' p.Results.inFileSuffix]));

% If we found at least one pupil data file, then proceed.
if ~isempty(fileListStruct)
    
    for ii = 1:length(fileListStruct)
        
        processThisOneFlag = true;
        % Check if we have a list of subjects to process
        if ~isempty(p.Results.subjectsToProcess)
            % Check if this subject is on the list
            if ~contains(fileListStruct(ii).name,p.Results.subjectsToProcess)
                processThisOneFlag = false;
            end
        end
        
        if processThisOneFlag
            % Convert this XML file to a volume
            mask = xml2volmask(fileListStruct(ii).name);
            
            % Create an output file path
            outFileRoot = extractBefore(fileListStruct(ii).name,p.Results.inFileSuffix);
            outputFile = [outFileRoot p.Results.outFileSuffix];
            
            % Save the file
            save(outputFile,'mask');
            
            % Report the conversion
            fprintf([outputFile '\n']);
        end
    end
end

end % main



function varargout = subdir(varargin)
%SUBDIR Performs a recursive file search
%
% subdir
% subdir(name)
% files = subdir(...)
%
% This function performs a recursive file search.  The input and output
% format is identical to the dir function.
%
% Input variables:
%
%   name:   pathname or filename for search, can be absolute or relative
%           and wildcards (*) are allowed.  If ommitted, the files in the
%           current working directory and its child folders are returned
%
% Output variables:
%
%   files:  m x 1 structure with the following fields:
%           name:   full filename
%           date:   modification date timestamp
%           bytes:  number of bytes allocated to the file
%           isdir:  1 if name is a directory; 0 if no
%
% Example:
%
%   >> a = subdir(fullfile(matlabroot, 'toolbox', 'matlab', '*.mat'))
%
%   a =
%
%   67x1 struct array with fields:
%       name
%       date
%       bytes
%       isdir
%
%   >> a(2)
%
%   ans =
%
%        name: '/Applications/MATLAB73/toolbox/matlab/audiovideo/chirp.mat'
%        date: '14-Mar-2004 07:31:48'
%       bytes: 25276
%       isdir: 0
%
% See also:
%
%   dir

% Copyright 2006 Kelly Kearney


%---------------------------
% Get folder and filter
%---------------------------

narginchk(0,1);
nargoutchk(0,1);

if nargin == 0
    folder = pwd;
    filter = '*';
else
    [folder, name, ext] = fileparts(varargin{1});
    if isempty(folder)
        folder = pwd;
    end
    if isempty(ext)
        if isdir(fullfile(folder, name))
            folder = fullfile(folder, name);
            filter = '*';
        else
            filter = [name ext];
        end
    else
        filter = [name ext];
    end
    if ~isdir(folder)
        error('Folder (%s) not found', folder);
    end
end

%---------------------------
% Search all folders
%---------------------------

pathstr = genpath_local(folder);
pathfolders = regexp(pathstr, pathsep, 'split');  % Same as strsplit without the error checking
pathfolders = pathfolders(~cellfun('isempty', pathfolders));  % Remove any empty cells

Files = [];
pathandfilt = fullfile(pathfolders, filter);
for ifolder = 1:length(pathandfilt)
    NewFiles = dir(pathandfilt{ifolder});
    if ~isempty(NewFiles)
        fullnames = cellfun(@(a) fullfile(pathfolders{ifolder}, a), {NewFiles.name}, 'UniformOutput', false);
        [NewFiles.name] = deal(fullnames{:});
        Files = [Files; NewFiles];
    end
end

%---------------------------
% Prune . and ..
%---------------------------

if ~isempty(Files)
    [~, ~, tail] = cellfun(@fileparts, {Files(:).name}, 'UniformOutput', false);
    dottest = cellfun(@(x) isempty(regexp(x, '\.+(\w+$)', 'once')), tail);
    Files(dottest & [Files(:).isdir]) = [];
end

%---------------------------
% Output
%---------------------------

if nargout == 0
    if ~isempty(Files)
        fprintf('\n');
        fprintf('%s\n', Files.name);
        fprintf('\n');
    end
elseif nargout == 1
    varargout{1} = Files;
end

end % subdir

function [p] = genpath_local(d)
% Modified genpath that doesn't ignore:
%     - Folders named 'private'
%     - MATLAB class folders (folder name starts with '@')
%     - MATLAB package folders (folder name starts with '+')

files = dir(d);
if isempty(files)
    return
end
p = '';  % Initialize output

% Add d to the path even if it is empty.
p = [p d pathsep];

% Set logical vector for subdirectory entries in d
isdir = logical(cat(1,files.isdir));
dirs = files(isdir);  % Select only directory entries from the current listing

for i=1:length(dirs)
    dirname = dirs(i).name;
    if    ~strcmp( dirname,'.') && ~strcmp( dirname,'..')
        p = [p genpath(fullfile(d,dirname))];  % Recursive calling of this function.
    end
end

end % genpath_local
