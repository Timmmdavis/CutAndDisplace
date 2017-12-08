function [varargout]=ReshapeData2d( rows,cols, varargin )
% ReshapeData2d: Presuming all the inputs are vectors or 2d matrices of the
%                   same size we can reshape them to the number of input
%                   rows and cols using this func. Created so you don't end
%                   up with multiple lines rehaping arrays such as tensors.
%               
% usage #1:
% [varargout]=ReshapeData2d( rows,cols, varargin )
%
% Arguments: (input)
% rows              - The number of rows the output matrix will have.
%
% cols              - The number of cols the output matrix will have.
%
% varargin          - Vectors/2d matricies to be reshaped.
%
% Arguments: (output)
% varargout         - As the varargin but reshaped to the input dimensions.
% 
%
% Example usage:
%
% %6*6 matrix 'X'
% X = (1:6)'*(2:7);
% Y = X*2; 
% Z=Y(:); 
% rows=2; cols=18;
% [X,Y,Z]=ReshapeData2d( rows,cols, X,Y,Z );
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%Preallocate cell array 
varargout= cell(1,nargout); 
%Start loop
for i=1:nargin-2
    %Getting inputs
    Data=varargin{i};
    
    %Checking the arrays are the right sizes
    DivTest=~mod(numel(Data),rows);
    if DivTest==0
        CurrentInputNm=inputname(i+2);
        %Reporting error for current input
        Str = strcat('Input array "',CurrentInputNm,'" cannot be reshaped to the sizes specified');
        error(Str)
    end
    
    %The actual reshape
    Data = reshape(Data,rows,cols); 
    %Assigning to outputs
    varargout{i}=Data;
end

