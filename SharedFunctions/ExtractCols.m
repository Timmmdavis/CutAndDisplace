function [ varargout ] = ExtractCols( Data )
% ExtractCols: Put in the input array and the amount of columns you want to
%                   take as the output args. The function does the rest. If
%                   you would like the first 3 columns of the input data
%                   specify 3 output arguments. 
%   
% usage #1:
% [ varargout ] = ExtractCols( Data )
%
% Arguments: (input)
% Data              - Large array that has columns. 
%
% Arguments: (output)
% varargout       	- Each column one by one extracted as a vector. 
%
% Example usage 1:
%
% %Grabbing each column from a matrix 'a':
% n=1000; 
% a=[rand(n,1),zeros(n,1),ones(n,1)];
% %Grab each column of a. 
% [ c1,c2,c3 ] = ExtractCols( a )
% %If you just wanted the 2nd row of a
% [ ~,c2,~ ] = ExtractArraysFromVector( a )
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


%Preallocate cell array
varargout= cell(1,nargout); 

%Start loop extracting the data:
for i=1:nargout
    %Getting the current first column
    varargout{i}=Data(:,1);
    %Remove the first column (this is more memory efficient) 
    Data=Data(:,2:end); 
end

end

