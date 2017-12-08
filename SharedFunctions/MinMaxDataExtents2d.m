function [maxgriX,mingriX,maxgriY,mingriY,size]=MinMaxDataExtents2d(points,cells,padding)
% MinMaxDataExtents2d: Creates limits of data for use in the function
%                   meshgrid. The limits are based on the xy limits of the
%                   input points and the resultant grid will be square,
%                   uniformly spaced in XY and have the number of cells
%                   defined by the input argument. Additional padding past
%                   the data limits can also be supplied.
%                   
%               
% usage #1:
% [maxgriX,mingriX,maxgriY,mingriY,size]=MinMaxDataExtents2d(points,cells,padding)
%
% Arguments: (input)
% Points           - XY locations of points in space. n*2 vector. 
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
% maxgriY,mingriY  - Data limits for use in the meshgrid function. 
%
% size             - Step size for use in the meshgrid function. 
%
% Example usage:
% 
% % Create a diagonal line that extends from -15 to 15 in X & Y.
% Points=[(-5:1:4)',(-4:1:5)',(-5:1:4)',(-4:1:5)']*3;
% % Plot this
% line([Points(:,1)';Points(:,2)'],[Points(:,3)';Points(:,4)'],'color','r')
% % Create a vector that is just all the points in X Y. 
% PointsXY=[[Points(:,1);Points(:,2)],[Points(:,3);Points(:,4)]]
% % Define number of desired cells and padding
% cells=10;
% padding=0;
% [maxgriX,mingriX,maxgriY,mingriY,sz]=MinMaxDataExtents2d(PointsXY,cells,padding);
% % Use results in MATLABs function meshgrid
% [X,Y] = meshgrid(mingriX:sz:maxgriX,mingriY:sz:maxgriY);
% % Draw this
% hold on
% scatter(X(:),Y(:));
% %Change padding
% padding=5;
% [maxgriX,mingriX,maxgriY,mingriY,sz]=MinMaxDataExtents2d(PointsXY,cells,padding);
% % Use results in MATLABs function meshgrid
% [X,Y] = meshgrid(mingriX:sz:maxgriX,mingriY:sz:maxgriY);
% % Draw this
% hold on
% scatter(X(:),Y(:),'r.');
% 
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%Getting the X extents of the data.
[maxx,minx,meanx]=Extents(points(:,1));
widthx=maxx-minx;
%Getting the Y extents of the data.
[mayy,miny,meany]=Extents(points(:,2));
widthy=mayy-miny;

%Now checking which is larger and building limits based on the data with
%the larger width.
if widthx>widthy
    size=widthx/cells;
    [maxgriX,mingriX,maxgriY,mingriY]=ExtentsFun(meanx,meany,widthx,size,padding);
else
    size=widthy/cells;
    [maxgriX,mingriX,maxgriY,mingriY]=ExtentsFun(meanx,meany,widthy,size,padding);
end

function [maxa,mina,meana]=Extents(a)
%Finds the limits and centre of the data.
maxa=max(a);
mina=min(a);
meana=mean(a);

function [maxgriX,mingriX,maxgriY,mingriY]=ExtentsFun(meanx,meany,width,size,padding)
%Calculates the grid limits for the data (padding included). 

% % % %For padding as a distance
% maxgriX=meanx+(width/2)+padding;
% mingriX=meanx-(width/2)-padding;
% maxgriY=meany+(width/2)+padding;
% mingriY=meany-(width/2)-padding;

% %For padding as no of additional cells:
maxgriX=meanx+(width/2)+(padding*size);
mingriX=meanx-(width/2)-(padding*size);
maxgriY=meany+(width/2)+(padding*size);
mingriY=meany-(width/2)-(padding*size);


