function [varargout]=ReshapeData3d( rows,cols,OthAx, varargin )
% ReshapeData3d: Presuming all the inputs are vectors or matrices of the
%                   same size we can reshape them to the number of input
%                   rows and cols and other size using this func. Created
%                   so you don't end up with multiple lines rehaping arrays
%                   such as tensors.
%               
% usage #1:
% [varargout]=ReshapeData2d( rows,cols, varargin )
%
% Arguments: (input)
% rows              - The number of rows the output matrix will have.
%
% cols              - The number of cols the output matrix will have.
%
% OthAx             - 3rd matrix dimension the output matrix will have.
%
% varargin          - Vectors/2d matricies to be reshaped.
%
% Arguments: (output)
% varargout         - As the varargin but reshaped to the input dimensions.
% 
%
% Example usage:
%
% X = 1:100;
% Y = X*2; 
% rows=2; cols=2; outplane=25;
% [X,Y]=ReshapeData3d( rows,cols,outplane, X,Y );
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


%Preallocate cell array 
varargout= cell(1,nargout); 
for i=1:nargin-3
    %Getting inputs
    Data=varargin{i};
    
    %Checking the arrays are the right sizes
    DivTest=numel(Data)==(rows*cols*OthAx);
    if DivTest==0
        CurrentInputNm=inputname(i+3);
        %Reporting error for current input
        Str = strcat('Input array "',CurrentInputNm,'" cannot be reshaped to the sizes specified');
        error(Str)
    end
    
    %The actual reshape
    Data = reshape(Data,rows,cols,OthAx); 
    %Assigning to outputs    
    varargout{i}=Data;
end

