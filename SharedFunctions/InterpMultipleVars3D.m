function [ varargout ] = InterpMultipleVars3D(X,Y,Z,Xq,Yq,Zq,varargin )
% InterpMultipleVars3D: Same as using MATLAB's interp3 but you can call
%                   multiple arguments for 'V'. Useful when interpolating
%                   tensors.
%
% usage #1:
% [ varargout ] = InterpMultipleVars3D(X,Y,Z,Xq,Yq,Zq,varargin )
%
% Arguments: (input)
% X,Y,Z             - Grid vectors XYZ, These vectors define the points
%                     associated with values in V (varargin).
%
% Xq,Yq,Zq          - Sample value location you plan of finding V at. 
%
% varargin          - The value you are interpolating. Must match the size
%                     of X,Y,Z.
%
% Arguments: (output)
% varargout         - Your value at points Xq,Yq,Zq.
%
% Example usage:
%
% [X,Y,Z]=meshgrid(-2:0.1:2,-2:0.1:2,-2:0.1:2);
% Density=5000;
% n=4;
% xmv=-2; ymv=-2;zmv=-2;
% Xq=(rand(1,Density)*n)'+xmv;
% Yq=(rand(1,Density)*n)'+ymv;
% Zq=(rand(1,Density)*n)'+zmv;
% % Interploating values of X, X^2 and X^3 at the random points (q). 
% [ Xq1,Xq2,Xq3 ] = InterpMultipleVars3D(X,Y,Z,Xq,Yq,Zq,X,X.^2,X.^3 );
% %Compare the two figs:
% %Figure of the interpolated values of X^2 at Xq.
% scatter3(Xq,Yq,Zq,[],Xq2);
% %Figure of actual values of Xq^2 at the random points.
% figure;
% scatter3(Xq,Yq,Zq,[],Xq.^2);
%
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


InputsSize=nargin-6;

for i=1:InputsSize
    varargin{i}= interp3(X,Y,Z,varargin{i},Xq,Yq,Zq);
end
varargout=varargin;


end

