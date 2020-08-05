function [Points,Triangles]=OFFReader(filename)

%Doing this as Octaves textscan is different to Matlabs
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0; %1 for octave, 0 for MATLAB
if  isOctave==1
disp('Can try octave but textscan is different')
elseif isOctave==0
end    

%% Initialize variables.
delimiter = ' ';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Format for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
textscan(fileID, '%[^\n\r]', startRow(1)-1, 'WhiteSpace', '', 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'EmptyValue', NaN, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    textscan(fileID, '%[^\n\r]', startRow(block)-1, 'WhiteSpace', '', 'ReturnOnError', false);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'EmptyValue', NaN, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Allocate imported array to column variable names
Col1 = dataArray{:, 1};
Col2 = dataArray{:, 2};
Col3 = dataArray{:, 3};
Col4 = dataArray{:, 4};

Collated=[Col1, Col2, Col3, Col4];

%Number of tris
n_points=Col1(1);
%Number of points
n_tris=Col2(1);
%Number of edges
n_edges=Col3(1);

Points=Collated(2:1+n_points,1:3);
Triangles=Collated(2+n_points:1+n_points+n_tris,2:4)+1;%ordering in .off starts at 0!

%Adding row with row numbers to the front of this: 
Sz=length(Points(:,1)); %getting size of rows
mono_inc=1:1:Sz; %monotomic increasing vec
Points=[mono_inc', Points];

end