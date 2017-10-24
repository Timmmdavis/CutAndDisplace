function [varargout]=ReshapeMultipleArrays(dimx,dimy,varargin)
%ReshapeMultipleArrays Pass as many arrays in as you want and reshape these
%to size dimx dimy. (only 2D reshaping) 

%   Copyright 2017, Tim Davis, The University of Potsdam

%Ignoring dimx and dimy inputs
InputsSize=nargin-2; 

for i=1:InputsSize
    varargin{i}=reshape(varargin{i},dimx,dimy);
end
varargout=varargin;

end
