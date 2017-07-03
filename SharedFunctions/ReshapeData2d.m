function [varargout]=ReshapeData2d( rows,cols, varargin )
%ReshapeData Presuming all the inputs are vectors of the same size we can
%reshape them to the number of input rows and cols using this func.
%For 2d matrices. Not reshaping in multidimensions here

% rows: the number of rows we want to reshape too
% cols: the number of cols we want to reshape too
% varargin: vectors/data we want to reshape

%   Copyright 2017, Tim Davis, The University of Aberdeen

varargout= cell(1,nargout); %Preallocate cell array 
for i=1:nargin-2
    Data=varargin{i};
    Data = reshape(Data,rows,cols);  %Sz=1=(number of Rows down),2=(number of Cols across)
    varargout{i}=Data;
end

