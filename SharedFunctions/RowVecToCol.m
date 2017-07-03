function [ varargout ] = RowVecToCol( varargin )
%RowToCol Checks vector and changes to column if it is a row

%varargin = chuck in as many arrays as you would like 
%varargout = get the arrays back but as cols

%   Copyright 2017, Tim Davis, The University of Aberdeen


InputsSize=nargin;

for i=1:InputsSize
    data=varargin{i}; %extract from inputs
    if isrow(data) %check its a row and make sure its a vect
    varargin{i}=data';%flip if its a row 
    end
end

varargout=varargin;  %add to the outputs

end

