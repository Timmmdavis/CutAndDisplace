function [ varargout ] = ExtractArraysFromVector( Data )
% ExtractArraysFromVector: Put in input vector and the amount of arrays you
%                   want as the outputs, the function does the rest. If you
%                   know a vector has 3 specific parts just put 3 output args.
%   
% usage #1:
% [ varargout ] = ExtractArraysFromVector( Data )
%
% Arguments: (input)
% Data              - Large vector that will be diced up into the output
%                     arrays.
%
% Arguments: (output)
% varargout       	- The input vector seperated into the number of output
%                     arguments 
%
% Example usage 1:
%
% %Grabbing bits of a col vector:
% n=1000; 
% a=[rand(n,1);zeros(n,1);ones(n,1)];
% %Splits 'a' into 3 bits
% [ top,mid,base ] = ExtractArraysFromVector( a )
% %If you just wanted the middle part of 'a'
% [ ~,mid,~ ] = ExtractArraysFromVector( a )
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


if ~isvector(Data)
    error('This function requires vector data as the input argument')
end

a=numel(Data);
b=nargout;

%Divisibility Check
Div=~mod(a,b); %1 if divisible
if ~Div
    error('The input vector is not divisible by the number of output arguments')
end

%Size of the vectors
size=a/b; 

%Preallocate cell array
varargout= cell(1,nargout); 

%Get the first output arg
varargout{1}=Data(1:size); 

    for i=2:nargout
        ii=i-1; 
        varargout{i}=Data((size*ii)+1:(size*i));
    end

end

