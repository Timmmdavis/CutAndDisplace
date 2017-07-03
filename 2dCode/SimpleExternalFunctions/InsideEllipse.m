function [ External ] = InsideEllipse(x, y, rx, ry, tol_dist)
%InsideEllipse finds if XY points are inside an ellipse with defined
%semi-axes
% X and Y are the location of the point. 
% rx and ry are semi-major axis rx, semi-minor axis ry
% both aligned with the Cartesian plane.

%   Copyright 2017, Tim Davis, The University of Aberdeen

%Getting original size
Xorig=size(x);
%Turning into col vecs
x=x(:);y=y(:);

sz=numel(x);
External = zeros(sz,1) ; 

for i=1:sz
    if ((x(i)^2)/rx^2+(y(i)^2)/ry^2) > 1.0+tol_dist 
    External(i) = true ; 
    else 
    %Do nothing
    end  
end
External=reshape(External,Xorig);


end

