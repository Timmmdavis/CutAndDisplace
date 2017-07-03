function [ varargout ] = ScaleMatrices( scl,varargin )
%Call to multiply multiple input matrices by the values in scl

%scl = your scale factor
%varargin = your arrays you want to scale, as many as you want

%example:
%[ A,B,C,D,E ] = ScaleMatrices( 5,A,B,C,D,E )

%   Copyright 2017, Tim Davis, The University of Aberdeen

InputsSize=nargin-1;

varargout= cell(1,nargout); %Preallocate cell array

for i=1:InputsSize

    data=varargin{i};
    data=data*scl;
    varargout{i}=data;

end


end

