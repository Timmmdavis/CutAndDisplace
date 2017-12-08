function [ X,Y ] = CreateCircleXY( Sampling,Radius,Closed )
% CreateCircleXY: Creates a circle of points in XY 
% 				around the edge of a circle with a defined radius,
% 				resultant number of points is defined by the 
%			 	sampling argument.
%
%               
% usage #1:
% [ X,Y ] = CreateCircleXY( Sampling,Radius )
%
% Arguments: (input)
% Sampling         - The resultant number of XY points 
%
% Radius           - The radius of the circle the points lie on. 
%
% Closed 	       - If the end of the circle is closed 
%				     (i.e. a duplicate point at its top).
%
% Arguments: (output)
% X,Y              - The resultant XY locations of points around the cicle. 
%
% Example usage 1:
% 
% Sampling=20;
% Radius=1;
% [ X,Y ] = CreateCircleXY( Sampling,Radius )
% plot(X,Y); hold on; scatter(X,Y); text(0,1.1,'Closed Circle, 20pnts')
% Closed=0; 
% [ X,Y ] = CreateCircleXY( Sampling,Radius,Closed )
% plot(X+2.1,Y); scatter(X+2.1,Y); text(2.1,1.1,'Open Circle, 20pnts')
% Closed=0; 
% [ X,Y ] = CreateCircleXY( Sampling-1,Radius,Closed )
% plot(X-2.1,Y); scatter(X-2.1,Y); text(-2.1,1.1,'Open Circle, 19pnts')
% axis('equal')
% 
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

if nargin==2
	Closed=1;
end

%Makes sure its the correct amount of points
if Closed==1
    Sampling=Sampling-1;
end

%Init vars
X=zeros(Sampling,1);
Y=X;



k = 0;
for theta = 0:pi/(Sampling/2):2*pi
	k = k+1;
	X(k) = Radius*sin(theta);
	Y(k) = Radius*cos(theta);
end

%Not closed, no duplicate point at top of circle
if Closed==0
	X(end)=[];
	Y(end)=[];
end

end

