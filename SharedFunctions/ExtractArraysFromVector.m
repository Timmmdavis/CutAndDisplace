function [ varargout ] = ExtractArraysFromVector( Data )
%ExtractArraysFromVector put in input array and the amount of arrays you
%want at the end, the function does the rest. If you know a vector has 3
%specific parts just put these as the output args. 

%   Copyright 2017, Tim Davis, The University of Aberdeen

if ~isvector(Data)
    error('This needs to be a vector for this func')
end

a=numel(Data);
b=nargout;

%Divisibility Check
Div=~mod(a,b); %1 if divisible
if ~Div
    error('Vector is not divisible by number of outputs')
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

