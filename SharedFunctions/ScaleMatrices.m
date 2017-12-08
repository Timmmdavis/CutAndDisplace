function [ varargout ] = ScaleMatrices( scl,varargin )
% ScaleMatrices: Call to multiply multiple input matrices by a scalar
%               
% usage #1:
% [ varargout ] = ScaleMatrices( scl,varargin )
%
% Arguments: (input)
% scl            - Scale factor 
%
% varargin       - chuck in as many arrays as you would like. 
%
% Arguments: (output)
% varargout      - get the arrays back but scaled
% 
% Example usage (1):
%
% A=ones(1,50);
% B=ones(50,1);
% Scale=700;
% [ A,B ] = ScaleMatrices( Scale,A,B );
% 
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%Get size of the input arguments. 
InputsSize=nargin-1;

for i=1:InputsSize
    %extract from inputs
    data=varargin{i};
    %Scale data
    data=data*scl;
    %Assign back to its input argument
    varargin{i}=data;

end

%add results to the outputs
varargout=varargin;  

end

