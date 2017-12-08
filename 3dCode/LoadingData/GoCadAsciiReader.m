function [ Points,Triangles ] = GoCadAsciiReader( Fid )
% GoCadAsciiReader: Loads Gocad Ascii files and extracts the points
%                   (triangle corners) and triangles (linked to points) of
%                   a triangulated mesh surface.
%
%                   This file makes the assumption your fault file in Fid
%                   is on one of MATLAB's current paths. 
%                   
%               
% usage #1:
% [ Points,Triangles ] = GoCadAsciiReader( Fid,pathstring )
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
%  string='SphereUniformDistribution5120FacesSides.ts';
%  [ Points,Triangles ] = GoCadAsciiReader( string );
%  PlotSlipDistribution3d(Triangles,Points,[],ones(size(Triangles(:,1))))
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

CurDir=pwd;
%'textscan' doesnt work with '.ts' files so we change these to '.dat'
[~, name, ~] =fileparts(Fid);
%Gets the location of the fault data (assuming its on the loaded paths) 
Loc=which(Fid); 
if isempty(Loc)
    error('Are you sure this file exists in your path?')
end
if ~isunix %Need to add, otherwise fails in linux+mac
    parts = strsplit(Loc, '\');
    DirPart = parts(1,1:end-1); %Removing the file name
    Loc = strjoin(DirPart,'\');
else %linux/mac has slashes the other way
    parts = strsplit(Loc, '/');
    DirPart = parts(1,1:end-1); %Removing the file name
    Loc = strjoin(DirPart,'/');
end

cd(Loc)
movefile(Fid, [name '.dat'])
Fid=[name '.dat'];

%% Format string for each line of text:
%   column1: text (%s)
%	column2: text (%s)
%   column3: text (%s)
%	column4: text (%s)
%   column5: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(Fid,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
delimiter = ' ';
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%Changing back to .ts
[~, name, ext] =fileparts(Fid);clear pathstr
tf = strcmp(ext,'.dat'); %Checking if the input string is a .ts, we need to rename to .dat for textscan to work
if tf==1
    %Only do if extension is .ts
    cd(Loc)
    %Renaming file back to .ts
    movefile(Fid, [name '.ts'])
    cd(CurDir) % back to this functions position
end

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
VarName1 = dataArray{:, 1};
VarName2 = dataArray{:, 2};
VarName3 = dataArray{:, 3};
VarName4 = dataArray{:, 4};
VarName5 = dataArray{:, 5};

%%
NAME = contains(VarName1,'PROJECTION');
Num = sum(NAME);
if Num>1
    error('More than one Gocad ascii surface in the file, this function can only deal with 1')
end    

%Flag if row is triangle or vertex, converting to logical so this can be
%used to extract rows
%Finding strings
VRTXflag = contains(VarName1,'VRTX');
TRGLflag = contains(VarName1,'TRGL');  
ATOMflag = contains(VarName1,'ATOM');  

%Converting to numeric before extracting data.
% VarName1 = str2double(VarName1);
VarName2 = str2double(VarName2);
VarName3 = str2double(VarName3);
VarName4 = str2double(VarName4);
VarName5 = str2double(VarName5);

%Extracting rows with points
Points=   [VarName2(VRTXflag),VarName3(VRTXflag),VarName4(VRTXflag),VarName5(VRTXflag)];
Triangles=[VarName2(TRGLflag),VarName3(TRGLflag),VarName4(TRGLflag)];
Atoms=    [VarName2(ATOMflag),VarName3(ATOMflag)]; %Duplicate verticies

%Now if atoms exist we need to plug these in to the points vector:
flag = any(Atoms);
if flag==1
    Tag=zeros(size(Atoms));
    Atoms=[Atoms,Tag]; %Adding rows to be filled
    %Looping and filling in correct XYZ for the atom numbers
    for i=1:numel(Atoms(:,1))
        %Find where points row 1 matches atom num
        a=Points(:,1)==Atoms(i,2);
        %Push in correct row
        Atoms(i,2:4)=Points(a,2:4);   
    end
    %Now appending atoms onto points and sorting this by row 1
    Points=[Points;Atoms];
    Points=sortrows(Points);
end    


%Fixing a bug, this codebase relies on the fact the cols are the same as
%numbers the triangles are pointing to in the first row of 'Points'.
%Some program exports do not get this right. 
for i=1:numel(Points(:,1))
    if Points(i,1) ~= i          %if the row is not equal to the first column on that row. 
        %any numbers correponding to this - 'Points(i,1)' need to be changed
        %to this 'i' in triangles
        NumToChange=Points(i,1);
        ToChangeTo=i;
        Bad = Triangles == NumToChange;   
        Triangles(Bad)=ToChangeTo;
        Points(i,1)=ToChangeTo;
    end
end   

cd(CurDir)%getting back to your original folder

end



