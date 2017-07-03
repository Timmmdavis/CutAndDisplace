function [ varargout ] = RepeatQuadValues( varargin )
%RepeatQuadValues Repeats values for quadratic elements. 
%   Code was setup for linear els with one collation point. Now we have 2
%   additional points on each el. In the coeff matrices etc these 2 new rows are
%   above and below the original collation points row. We repeat matrices
%   such as normal directions using this func. 

%Input such as [1;2;3]
%gives an
%Output = [1;1;1;2;2;2;3;3;3]

InputsSize=nargin;

for i=1:InputsSize
    varargin{i}=[varargin{i}';varargin{i}';varargin{i}'];
    varargin{i}=varargin{i}(:);
end

varargout=varargin;
end

