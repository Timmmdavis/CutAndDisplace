function  DrawStressEllipsoidsPrincipal(StressTensor,StrainTensor,X,Y,Z,varargin)
% DrawStressEllipsoidsPrincipal: Draws the principal stress ellipsoids, 
%                   pretty intensive function so can be quite slow. A catch
%                   exists to stop users putting too much data into this
%                   and taking forever.
%
%                   Results show as red ellipsoids when the dilatation is
%                   positive and blue when negative.
%
%                   This also draws a triangulated surface if you supply
%                   this in the varargin. 
%               
% usage #1:
% DrawStressEllipsoidsPrincipal(StressTensor,StrainTensor,X,Y,Z,'Scale',2)
%
% Arguments: (input)
% StressTensor      - From the stress calculation 6*n vector, 
%                    [Sxx,Syy,Szz,Sxy,Sxz,Syz]
%
% StrainTensor      - From the strain calculation 6*n vector, 
%                    [Exx,Eyy,Ezz,Exy,Exz,Eyz]
%
% XYZ               - The XYZ locations of the points where the tensors are.
%
% Additional arguments: (input)
%
% Scale             - Property changing the length of the 
%                    resultant ellipsoid. (Single value)  
%
% Sample            - Properties changing the sampling of the 
%                    resultant ellipsoid. (Single value). Its the amount of
%                    points/faces that make up each ellipsoid in plot, the
%                    points are Sample^2. Reduce this number if figure is
%                    taking too long to draw.
%
% Points            - Columns 2 3 and 4 are the XYZ locations of one the
%                    corner points of a triangle. Column 1 is the index.
%                    Not needed unless you want to draw a surface too. 
%
% Triangles         -  Triangles is a list where each row contains 3 index
%                    locations in "Points" which contains the XYZ location
%                    of each corner of the triangle.
%                    Not needed unless you want to draw a surface too. 
%
% Arguments: (output)
% N/A - This just draws a figure. 
%
% Example usage:
%
% %Random Points on a sphere. 
% Theta = 2*pi*rand(1,250);
% Phi = asin(-1+2*rand(1,250));
% [X,Y,Z] = sph2cart(Theta',Phi',1);
% Sxx=ones(size(X)).*Theta';
% Syy=ones(size(X)).*-Theta';
% Szz=ones(size(X)).*Phi';
% Sxy=Syy;
% Sxz=Syy;
% Syz=Syy;
% [ Exx,Eyy,Ezz,Exy,Exz,Eyz ] = HookesLaw3dStress2Strain( Sxx,Syy,Szz,Sxy,Sxz,Syz,5,5 );
% StrainTensor=[Exx,Eyy,Ezz,Exy,Exz,Eyz]
% StressTensor=[Sxx,Syy,Szz,Sxy,Sxz,Syz];
% DrawStressEllipsoidsPrincipal(StressTensor,StrainTensor,X,Y,Z,'Scale',0.01)
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


%Checking if additional arguments have been input into the function:
[ varargin,Scale ]  = AdditionalArgsInVaragin( 'Scale',     varargin,1 );
[ varargin,Sample ] = AdditionalArgsInVaragin( 'Sample',    varargin,15 );
[ varargin,Points ] = AdditionalArgsInVaragin( 'Points',    varargin,[] );
[ ~,Triangles ]     = AdditionalArgsInVaragin( 'Triangles', varargin,[] );


%Loading external Octave packages, Normc etc are in here
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0; %1 for octave, 0 for MATLAB
if  isOctave==1
    pkg load general
    pkg load miscellaneous
    elseif isOctave==0
    %do nothing
end

if numel(X)>1000
    disp('Leaving function "DrawStressEllipsoidsPrincipal.m", too many input points')
    disp('It is recommended to use a few points at first to get the parameters right.')
    return
end

figure;

%Forcing XYZ to column vecs 
Xcv=X(:); %col vecs
Ycv=Y(:); %col vecs
Zcv=Z(:); %col vecs

%Extracting strain & stress
[Sxx,Syy,Szz,Sxy,Sxz,Syz ] = ExtractCols( StressTensor );
[Exx,Eyy,Ezz,Exy,Exz,Eyz ] = ExtractCols( StrainTensor );
%Calculating principal stresses
[S1,S2,S3,S1dir,S2dir,S3dir] = EigCalc3d(Sxx(:),Syy(:),Szz(:),Sxy(:),Sxz(:),Syz(:));
%Calculating principal strains to calculate dilatation. (if vol change is
%positive or negative for drawing
[E1,E2,E3] = EigCalc3d(Exx(:),Eyy(:),Ezz(:),Exy(:),Exz(:),Eyz(:));
Dilatation=E1+E2+E3; %(Change/Orig vol)
Expanded=Dilatation(:)>0; %flag, positive is expansion

%Init Dims Outside lp
dimx = Sample+1;
dimy = Sample+1;  

%Now drawing first subplot
%This is a figure showing vectors that are the principal directions.
%Coloured for S1S2S3 and sign, see fig title for col details 
for i=1:numel(S1)
    %Creating the ellipsoid
    [x,y,z] = ellipsoid(0,0,0,S1(i)*Scale,S2(i)*Scale,S3(i)*Scale,Sample); %Very low sampling
    
    %Rotating points to each direction cosine
    [xn,yn,zn] = RotateObject3dNewCoords(S1dir(i,:),S2dir(i,:),S3dir(i,:),x,y,z);

    %Reshape
    [xn,yn,zn]=ReshapeData2d( dimx,dimy, xn,yn,zn );
    
    
    %Moving the points to the correct location in space
    xn=xn+Xcv(i,:);
    yn=yn+Ycv(i,:);
    zn=zn+Zcv(i,:);

    [ C ] = RGB2Colour(0,191,255);

        %Drawing a single ellipsoid then setting the colour based on the total volume change sign at this
        %infinitesimal point. 
        if Expanded(i)==1
            hSurface=surf(xn, yn, zn);
            if  isOctave==1
                set(hSurface,'FaceColor','red','EdgeColor',([0.3;0.3;0.3]));%lighting not setup for surfaces in octave
            elseif isOctave==0
                set(hSurface,'FaceColor','red','FaceLighting','gouraud','EdgeColor',([0.3;0.3;0.3])); 
            end
        else
            hSurface=surf(xn, yn, zn);
            if  isOctave==1
                set(hSurface,'FaceColor',C,'EdgeColor',([0.3;0.3;0.3]));%lighting not setup for surfaces in octave
            elseif isOctave==0
                set(hSurface,'FaceColor',C,'FaceLighting','gouraud','EdgeColor',([0.3;0.3;0.3])); 
            end
        end    
    hold on
end
xlabel('x'); ylabel('y');zlabel('z');axis equal ;

%Only setting tranparency if MATLAB
if isOctave==0
    alpha(.8);%transparency value
end
title({'\fontsize{14}Stress ellipsoids','\fontsize{8}Blue, net compressional volume change, Red, Net extensional volume change'})

if ~isempty(Triangles)
    hold on
    trisurf(Triangles,Points(:,2),Points(:,3),Points(:,4),'FaceAlpha',(.2),'FaceColor', [0.5 0 0.9 ]);
    hold off
end

end

