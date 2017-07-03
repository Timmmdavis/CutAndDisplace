%   Copyright 2017, Tim Davis, The University of Aberdeen

function [maxgriX,mingriX,maxgriY,mingriY,size]=MinMaxDataExtents(points,cells,padding)

[maxx,minx,meanx]=Extents([points(:,1);points(:,2)]);
extentsx=maxx-minx;
[mayy,miny,meany]=Extents([points(:,3);points(:,4)]);
extentsy=mayy-miny;
if extentsx>extentsy
    size=extentsx/cells;
    [maxgriX,mingriX,maxgriY,mingriY]=ExtentsFun(meanx,meany,extentsx,size,padding);
else
    size=extentsy/cells;
    [maxgriX,mingriX,maxgriY,mingriY]=ExtentsFun(meanx,meany,extentsy,size,padding);
end

function [maxa,mina,meana]=Extents(a)
maxa=max(a);
mina=min(a);
meana=mean(a);

function [maxgriX,mingriX,maxgriY,mingriY]=ExtentsFun(meanx,meany,extents,size,padding)
maxgriX=meanx+(extents/2)+(padding*size);
mingriX=meanx-(extents/2)-(padding*size);
maxgriY=meany+(extents/2)+(padding*size);
mingriY=meany-(extents/2)-(padding*size);



