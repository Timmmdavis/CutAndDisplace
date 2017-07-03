function [ Points,Triangles ] = GoCadAsciiReader( Fid,pathstring )
%GoCadAsciiReaderg
%Loads Gocad Ascii files and extracts the points and triangles this
%software uses. 
%Fid is a string containing the file ID

%   Copyright 2017, Tim Davis, The University of Aberdeen
%Doing this as Octaves Textscan is hard to get right, using an old function
%I wrote that is not as solid
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0; %1 for octave, 0 for MATLAB
if  isOctave==1;
[ Points,Triangles ] = GoCadAsciiReaderOctave( Fid );
elseif isOctave==0;


%TEXTSCAN DOES NOT WORK WITH .TS, changing to dat
[pathstr, name, ext] =fileparts(Fid);clear pathstr
tf = strcmp(ext,'.ts'); %Checking if the input string is a .ts, we need to rename to .dat for textscan to work
if tf==1;
%Only do if extension is .ts
%===== Getting to Fault Data directory ==========================
AdRs0=mfilename('fullpath'); %Directory we are running from

if ~isunix %Need to add, otherwise fails in linux+mac
parts = strsplit(AdRs0, '\');
DirPart = parts(1,1:end-1); %Removing the file name
AdRs0 = strjoin(DirPart,'\');
else %linux has slashes the other way
parts = strsplit(AdRs0, '/');
DirPart = parts(1,1:end-1); %Removing the file name
AdRs0 = strjoin(DirPart,'/');
end
cd(AdRs0) %Making sure we are in the place we are running from.
cd('..') % go up one level, we are in folder 1optB
cd FaultData %go into the fault data directoy
%Renaming file
movefile(Fid, [name '.dat'])
Fid=[name '.dat'];
%cd(AdRs0) % back to your original position WRONG
end

delimiter = ' ';

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
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%Changing back to .ts
[pathstr, name, ext] =fileparts(Fid);clear pathstr
tf = strcmp(ext,'.dat'); %Checking if the input string is a .ts, we need to rename to .dat for textscan to work
if tf==1;
%Only do if extension is .ts
cd(AdRs0) %Making sure we are in the place we are running from.
cd('..') % go up one level, we are in folder 1optB
cd FaultData %go into the fault data directoy
%Renaming file back to .ts
movefile(Fid, [name '.ts'])
cd(AdRs0) % back to this functions position
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

NAME = strfind(VarName1, 'PROJECTION');
NAME = ~cellfun(@isempty,NAME);
Num = sum(NAME);
if Num>1;
    error('More than one Gocad ascii surface in the file, this function can only deal with 1')
end    

%Flag if row is triangle or vertex, converting to logical so this can be
%used to extract rows
%Finding strings
VRTXflag = strfind(VarName1, 'VRTX');
TRGLflag = strfind(VarName1, 'TRGL');
ATOMflag = strfind(VarName1, 'ATOM');
%Replacing empty cells with zeros, this converts the cell data to a logical
VRTXflag = ~cellfun(@isempty,VRTXflag);
TRGLflag = ~cellfun(@isempty,TRGLflag);
ATOMflag = ~cellfun(@isempty,ATOMflag);


%Converting to numeric before extracting data.
VarName1 = str2double(VarName1);
VarName2 = str2double(VarName2);
VarName3 = str2double(VarName3);
VarName4 = str2double(VarName4);
VarName5 = str2double(VarName5);

%Extracting rows with points
Points=[VarName2(VRTXflag),VarName3(VRTXflag),VarName4(VRTXflag),VarName5(VRTXflag)];
Triangles=[VarName2(TRGLflag),VarName3(TRGLflag),VarName4(TRGLflag)];
Atoms=[VarName2(ATOMflag),VarName3(ATOMflag)]; %Duplicate verticies

%Now if atoms exist we need to plug these in to the points vector:
flag = any(Atoms);
if flag==1
    Tag=zeros(size(Atoms));
    Atoms=[Atoms,Tag]; %Adding rows to be filled
    %Looping and filling in correct XYZ for the atom numbers
    for i=1:numel(Atoms(:,1));
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
%Some program exports don’t get this right. 
for i=1:numel(Points(:,1));
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

end
cd(pathstring)%getting back to your original folder
end

