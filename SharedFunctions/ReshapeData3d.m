function [varargout]=ReshapeData3d( rows,cols,OthAx, varargin )
%ReshapeData Presuming all the inputs are vectors of the same size we can
%reshape them to the number of input rows and cols using this func

% rows: the number of rows we want to reshape too
% cols: the number of cols we want to reshape too
% OthAx: 3rd matrix dimension. 'other axis'
% varargin: vectors/data we want to reshape

%   Copyright 2017, Tim Davis, The University of Aberdeen

varargout= cell(1,nargout); %Preallocate cell array 
for i=1:nargin-3
    Data=varargin{i};
    Data = reshape(Data,rows,cols,OthAx);  %Sz=1=(number of Rows down),2=(number of Cols across)
    varargout{i}=Data;
end

