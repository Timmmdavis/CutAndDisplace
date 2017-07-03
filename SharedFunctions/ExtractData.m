function [ varargout ] = ExtractData( Flag,type,varargin )
%ExtractData just grabs columns or rows specified by a flag. Saves space and makes
%more important functions clear

%   Copyright 2017, Tim Davis, The University of Aberdeen


%type is a number if:
%1=Take rows specified by flag
%2=Remove rows specified by flag
%1=Take Cols specified by flag
%2=Remove Cols specified by flag

InputsSize=nargin-2;

for i=1:InputsSize
    if type==1
    data=varargin{i};
    data=data(Flag,:);
    varargin{i}=data;
    elseif type==2
    data=varargin{i};
    data=data(~Flag,:);
    varargin{i}=data;    
    elseif type==3
    data=varargin{i};
    data=data(:,Flag);
    varargin{i}=data;    
    elseif type==4
    data=varargin{i};
    data=data(:,~Flag);
    varargin{i}=data;
    end
end

varargout=varargin;



end

