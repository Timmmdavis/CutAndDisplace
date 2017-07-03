function [varargout]=FlipColVecs(varargin)
%FlipColVects Pass as many as you want in and get the same amount out, just flipped. 
%They can be different sizes

%   Copyright 2017, Tim Davis, The University of Aberdeen


InputsSize=nargin;

for i=1:InputsSize
    varargin{i}=flipud(varargin{i});
end

varargout=varargin;

end

