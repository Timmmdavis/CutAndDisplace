%   Copyright 2017, Tim Davis, The University of Aberdeen
function [maxgriX,mingriX,maxgriY,mingriY,maxgriZ,mingriZ,size]=MinMaxDataExtents3d(points,cells,padding)

[maxx,minx,meanx]=Extents(points(:,2));
extentsx=maxx-minx;
[maxy,miny,meany]=Extents(points(:,3));
extentsy=maxy-miny;
[maxz,minz,meanz]=Extents(points(:,4));
extentsz=maxz-minz;

%placing in vector
arr = [extentsx,extentsy,extentsz];
[maximum,ind] = max(arr);


if ind==1
    size=extentsx/cells;
    [maxgriX,mingriX,maxgriY,mingriY,maxgriZ,mingriZ]=ExtentsFun(meanx,meany,meanz,extentsx,size,padding);
elseif ind==2
    size=extentsy/cells;
    [maxgriX,mingriX,maxgriY,mingriY,maxgriZ,mingriZ]=ExtentsFun(meanx,meany,meanz,extentsy,size,padding);
else  %must be z 
    size=extentsz/cells;
    [maxgriX,mingriX,maxgriY,mingriY,maxgriZ,mingriZ]=ExtentsFun(meanx,meany,meanz,extentsz,size,padding);
end


function [maxa,mina,meana]=Extents(a)
maxa=max(a);
mina=min(a);
meana=mean(a);

function [maxgriX,mingriX,maxgriY,mingriY,maxgriZ,mingriZ]=ExtentsFun(meanx,meany,meanz,extents,size,padding)
% maxgriX=meanx+(extents/2)+(padding*size);
% mingriX=meanx-(extents/2)-(padding*size);
% maxgriY=meany+(extents/2)+(padding*size);
% mingriY=meany-(extents/2)-(padding*size);
% maxgriZ=meanz+(extents/2)+(padding*size);
% mingriZ=meanz-(extents/2)-(padding*size);
maxgriX=meanx+(extents/2)+padding;
mingriX=meanx-(extents/2)-padding;
maxgriY=meany+(extents/2)+padding;
mingriY=meany-(extents/2)-padding;
maxgriZ=meanz+(extents/2)+padding;
mingriZ=meanz-(extents/2)-padding;


