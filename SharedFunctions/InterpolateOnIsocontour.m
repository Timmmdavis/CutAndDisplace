function [ DATA_SURF_X,DATA_SURF_Y,DATA_SURF_Z,varargout ] = InterpolateOnIsocontour( X,Y,Z,Data,Value,varargin )
%InterpolateOnIsocontour 3d function to interpolate data onto a 3d surface
%or isocontour
%   Input data such as tensors, interpolate multiple tensors then outside this func you can draw vector fields of these tensors
%   and have these sitting on the isocontour, For example stress directions on a principal stress
%   contour surface. 

%   Copyright 2017, Tim Davis, The University of Aberdeen

%XYZ values on a grid
%Data - gridded data that you find the isocontour of, for example Stress-S1
%Value - the value at which you want the isocontour
%varargin - the tensors that you want interpolated onto the contour

%Interpolated data - DATA_SURF_X,DATA_SURF_Y,DATA_SURF_Z,varargout
%The xyz and the new interpolated tensors/data at these points

%creating the isosurface
SURF = isosurface(X,Y,Z,Data,Value);
%interpolating tensors at this value
%first XYZ of the new data
DATA_SURF_X = interp3(X,Y,Z,X,SURF.vertices(:,1),SURF.vertices(:,2),SURF.vertices(:,3));
DATA_SURF_Y = interp3(X,Y,Z,Y,SURF.vertices(:,1),SURF.vertices(:,2),SURF.vertices(:,3));
DATA_SURF_Z = interp3(X,Y,Z,Z,SURF.vertices(:,1),SURF.vertices(:,2),SURF.vertices(:,3));


varargout= cell(1,nargout); %Preallocate cell array 
for i=1:nargin-5
    Data=varargin{i};
    %and the interpolated data
    DATA_SURF_Interp = interp3(X,Y,Z,Data,SURF.vertices(:,1),SURF.vertices(:,2),SURF.vertices(:,3));    
    varargout{i}=DATA_SURF_Interp;
end


end

