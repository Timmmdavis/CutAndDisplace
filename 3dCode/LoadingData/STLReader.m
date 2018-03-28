function [ Points,Triangles ] = STLReader( Fid)
% STLReader: Imports Mesh STL files. Loads non binary STL files and
%                   extracts the points and triangles this software uses.
%                   Fid is a string containing the file ID. If you are
%                   having issues import the mesh into 'meshlab', export as
%                   a .stl and make sure in additional parameters binary
%                   encoding is turned OFF.
%                   
%               
% usage #1:
% [ Points,Triangles ] = STLReader( Fid )
%
% Arguments: (input)
% Fid               - Fid is a string containing the file ID.
%
% Arguments: (output)
% Points            - Columns 2 3 and 4 are the XYZ locations of one the
%                    corner points of a triangle. Column 1 is the index.
%                    Not needed unless you want to draw a surface too. 
%
% Triangles         -  Triangles is a list where each row contains 3 index
%                    locations in "Points" which contains the XYZ location
%                    of each corner of the triangle.
%                    Not needed unless you want to draw a surface too. 
%
%
% Example usage:
%   string='Sphere5120.stl';
%  [ Points,Triangles ] = STLReader( string );
%  PlotSlipDistribution3d(Triangles,Points,[],ones(size(Triangles(:,1))))
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


%Doing this as Octaves textscan is different to Matlabs
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0; %1 for octave, 0 for MATLAB
if  isOctave==1
disp('Can try octave but textscan is different')
elseif isOctave==0
end    

%% Initialize variables.
delimiter = ' ';
startRow = 2;

%% Read columns of data as strings:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%[^\n\r]';

%% Open the text file.
%fileID = fopen(filename,'r');
fileID = fopen(Fid,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Allocate imported array to column variable names
VarName1 = dataArray{:, 1};
VarName2 = dataArray{:, 2};
VarName3 = dataArray{:, 3};
VarName4 = dataArray{:, 4};

%Flag every row that is a vertex, in stls these come in groups of threes. 
if verLessThan('matlab', '9.1') %(below v2016b where 'contains' was introduced)
    [VRTXflag,~] = find(~cellfun(@isempty,strfind(VarName1,'vertex')));
else
    VRTXflag = contains(VarName1, 'vertex');
end
%Extracting rows with points
Points=[VarName2(VRTXflag),VarName3(VRTXflag),VarName4(VRTXflag)];
Points = cellfun(@str2double, Points); %converting to double

%Adding row with row numbers to the front of this: 
Sz=numel(Points(:,1)); %getting size of rows
mono_inc=(1:1:Sz)'; %monotomic increasing vec
Points=[mono_inc,Points];

%Creating triangles pointer
Triangles=1:numel(Points(:,1));
Sz2=numel(Triangles)/3;
Triangles=reshape(Triangles.',[],Sz2)';
