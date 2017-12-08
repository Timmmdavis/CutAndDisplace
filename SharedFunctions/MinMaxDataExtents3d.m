function [maxgriX,mingriX,maxgriY,mingriY,maxgriZ,mingriZ,size]=MinMaxDataExtents3d(points,cells,padding)
% MinMaxDataExtents3d: Creates limits of data for use in the function
%                   meshgrid. The limits are based on the xyz limits of the
%                   input points and the resultant grid will be square,
%                   uniformly spaced in XYZ and have the number of cells
%                   defined by the input argument. Additional padding past
%                   the data limits can also be supplied.
%                   
%               
% usage #1:
% [maxgriX,mingriX,maxgriY,mingriY,maxgriZ,mingriZ,size]=MinMaxDataExtents3d(points,cells,padding)
%
% Arguments: (input)
% Points           - XYZ locations of points in space. n*3 vector. 
%
% cells            - This creates gridded data, this defines how many
%                   cells span across the data extents. 
%
% padding          - The number of extra cells you would like away from the
%                   data limits (not distance). Using this means there will
%                   be more cells that defined in the input args. 
%
% Arguments: (output)
% maxgriX,mingriX
% maxgriY,mingriY  
% maxgriZ,mingriZ  - Data limits for use in the meshgrid function. 
%
% size             - Step size for use in the meshgrid function. 
%
% Example usage:
% 
% %Create triangulated surface
% [x,y] = meshgrid(-2:.2:2);                                
% z = x .* exp(-x.^2 - y.^2);
% Triangles = delaunay(x(:),y(:));
% Points=[x(:),y(:),z(:)];
% trisurf(Triangles,Points(:,1),Points(:,2),Points(:,3));
% % Define number of desired cells and padding
% cells=10;
% padding=0;
% [maxgriX,mingriX,maxgriY,mingriY,maxgriZ,mingriZ,sz]=MinMaxDataExtents3d(Points,cells,padding);
% % Use results in MATLABs function meshgrid
% [X,Y,Z] = meshgrid(mingriX:sz:maxgriX,mingriY:sz:maxgriY,mingriZ:sz:maxgriZ); 
% hold on
% scatter3(X(:),Y(:),Z(:))
% % Change padding
% padding=2;
% [maxgriX,mingriX,maxgriY,mingriY,maxgriZ,mingriZ,sz]=MinMaxDataExtents3d(Points,cells,padding);
% % Use results in MATLABs function meshgrid
% [X,Y,Z] = meshgrid(mingriX:sz:maxgriX,mingriY:sz:maxgriY,mingriZ:sz:maxgriZ); 
% hold on
% scatter3(X(:),Y(:),Z(:),1,'r')
% 
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%Getting the X extents of the data.
[maxx,minx,meanx]=Extents(points(:,1));
widthx=maxx-minx;
%Getting the Y extents of the data.
[maxy,miny,meany]=Extents(points(:,2));
widthy=maxy-miny;
%Getting the Z extents of the data.
[maxz,minz,meanz]=Extents(points(:,3));
widthz=maxz-minz;

%Placing in vector
arr = [widthx,widthy,widthz];
%Finding index of max width
[~,ind] = max(arr);

%Now checking which is larger and building limits based on the data with
%the larger width.
if ind==1
    size=widthx/cells;
    [maxgriX,mingriX,maxgriY,mingriY,maxgriZ,mingriZ]=ExtentsFun(meanx,meany,meanz,widthx,size,padding);
elseif ind==2
    size=widthy/cells;
    [maxgriX,mingriX,maxgriY,mingriY,maxgriZ,mingriZ]=ExtentsFun(meanx,meany,meanz,widthy,size,padding);
else  %must be z 
    size=widthz/cells;
    [maxgriX,mingriX,maxgriY,mingriY,maxgriZ,mingriZ]=ExtentsFun(meanx,meany,meanz,widthz,size,padding);
end


function [maxa,mina,meana]=Extents(a)
%Finds the limits and centre of the data.
maxa=max(a);
mina=min(a);
meana=mean(a);

function [maxgriX,mingriX,maxgriY,mingriY,maxgriZ,mingriZ]=ExtentsFun(meanx,meany,meanz,width,size,padding)
%Calculates the grid limits for the data (padding included). 

% % %For padding as a distance
% maxgriX=meanx+(width/2)+padding;
% mingriX=meanx-(width/2)-padding;
% maxgriY=meany+(width/2)+padding;
% mingriY=meany-(width/2)-padding;
% maxgriZ=meanz+(width/2)+padding;
% mingriZ=meanz-(width/2)-padding;

%For padding as no of additional cells:
maxgriX=meanx+(width/2)+(padding*size);
mingriX=meanx-(width/2)-(padding*size);
maxgriY=meany+(width/2)+(padding*size);
mingriY=meany-(width/2)-(padding*size);
maxgriZ=meanz+(width/2)+(padding*size);
mingriZ=meanz-(width/2)-(padding*size);



