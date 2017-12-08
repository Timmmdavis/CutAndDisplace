function [ XCont,YCont,ZCont,varargout ] = InterpolateOnIsocontour( X,Y,Z,Data,Value,varargin )
% InterpolateOnIsocontour: function to interpolate gridded data onto a 3d
%                           isocontour. Can take multiple arguments (such
%                           as tensor fields).
%
% usage #1:
% [ XCont,YCont,ZCont,varargout ] = InterpolateOnIsocontour( X,Y,Z,Data,Value,varargin )
%
% Arguments: (input)
% X,Y,Z             - Grid vectors XYZ, These vectors define the points
%                     associated with values in V (varargin).
%
% Data              - The data you are finding an isocontour on. 
%
% Value             - The value that defines the isocontour's location. 
%
% varargin          - The value you are interpolating. Must match the size
%                     of X,Y,Z.
%
% Arguments: (output)
%  XCont,YCont,ZCont- The locations of the interpolated data on the surface
%                     in 3D space.
%
% varargin          - The interpolated values of the input args at the 
%                     locations defined by XCont,YCont,ZCont.
%
% Example usage:
%
% %Create a grid of points
% [X,Y,Z]=meshgrid(-2:0.2:2, -2:0.2:2, -2:0.2:2);
% dimx = length(X(:,:,1)); 
% dimy = length(X(:,1,:));
% dimz = length(X(1,:,:));
% [~,~,r] = cart2sph(X(:),Y(:),Z(:));
% %Create displacements that represent inflation from 0,0. (Big bang style).
% %Further from origin = larger displacement. 
% Ux=(r(:).^2).*X(:);
% Uy=(r(:).^2).*Y(:);
% Uz=(r(:).^2).*Z(:);
% %Reshape this
% [Ux,Uy,Uz]=ReshapeData3d(dimx,dimy,dimz,Ux,Uy,Uz);
% %Show we have just have expansion from 0,0 increasing away from 0. 
% quiver3(X,Y,Z,Ux,Uy,Uz)
% %Compute strain tensors for this data
% [exx,eyy,ezz,exy,exz,eyz]=FiniteStrainLagrangian3d(Ux,Uy,Uz,X,Y,Z);
% %Compute eigen values for this on the grid
% [E1,E2,E3,E1dir,E2dir,E3dir] = EigCalc3d(exx,eyy,ezz,exy,exz,eyz);
% %Reshape results
% [X,Y,Z,r,exx,eyy,ezz,exy,exz,eyz]=ReshapeData3d(dimx,dimy,dimz,X,Y,Z,r,exx,eyy,ezz,exy,exz,eyz);
% %Draw the 50th percentile of the radius value (should be a ball).
% r_50 = Percentile(r(:),50); 
% figure;
% P1 = patch(isosurface(X,Y,Z,r,r_50)) ; 
% set(P1,'FaceColor', RGB2Colour(255, 0, 234),'EdgeColor','none') 
% camlight ; lighting gouraud ;
% hold on
% %Calculate the strain tensors on the contour of r @ 50th percentile
% [ x,y,z,exx,eyy,ezz,exy,exz,eyz  ]...
%     = InterpolateOnIsocontour( X,Y,Z,r,r_50,exx,eyy,ezz,exy,exz,eyz  );
% %Compute eigen values for these interpolated values on the contour of r. 
% [~,~,~,e1dir,e2dir,e3dir] = EigCalc3d(exx,eyy,ezz,exy,exz,eyz);
% %Now draw the directions of the most extensional strain on this contour.
% %This msut be perpendicular to the surface.
% quiver3(x,y,z,e1dir(:,1),e1dir(:,2),e1dir(:,3),'b')
% quiver3(x,y,z,-e1dir(:,1),-e1dir(:,2),-e1dir(:,3),'b')
% axis('equal')
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


%Creating the isosurface at the value defined by 'Value'
SURF = isosurface(X,Y,Z,Data,Value);
%Interpolating tensors at this value
%First grab XYZ of the new data
Xq=SURF.vertices(:,1); Yq=SURF.vertices(:,2); Zq=SURF.vertices(:,3);
%Doing with external func
[ XCont,YCont,ZCont ] = InterpMultipleVars3D(X,Y,Z,Xq,Yq,Zq,X,Y,Z );

%Preallocate cell array 
varargout= cell(1,nargout); 

for i=1:numel(varargin)
    Data=varargin{i};
    %and the interpolated data
    DATA_SURF_Interp = interp3(X,Y,Z,Data,SURF.vertices(:,1),SURF.vertices(:,2),SURF.vertices(:,3));    
    varargout{i}=DATA_SURF_Interp;
end


end

