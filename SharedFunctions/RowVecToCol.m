function [ varargout ] = RowVecToCol( varargin )
% RowVecToCol: Checks vector and changes to column if it is a row. Works
%               with multiple arguments.
%               
% usage #1:
% [ varargout ] = RowVecToCol( varargin )
%
% Arguments: (input)
% varargin       - chuck in as many arrays as you would like. 
%
% Arguments: (output)
% varargout      - get the arrays back but as cols.
% 
% Example usage (1):
%
% A=ones(1,50);
% B=ones(50,1);
% [ A,B ] = RowVecToCol( A,B )
% 
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%Get size of the input arguments. 
InputsSize=nargin;

for i=1:InputsSize
    %extract from inputs
    data=varargin{i}; 
    %check its a row and make sure its a vect
    if isrow(data) 
        %flip if its a row 
        varargin{i}=data';
    end
end
%add results to the outputs
varargout=varargin;  

end

