function [ x,y,xe,ye,HalfLength,Beta,CosB,Points,NormAng ] = MidPoint_Orientation( XBEG,XEND,YBEG,YEND,NUM )
% Extracting MidPoints, size and normal orientation for each element
% The outward normal is measured from Xaxis counter clockwise to fit with
% the convention described in Pollard and Fletcher Fundamentals 2005
% Chapter 6.2 Page 215

%   Copyright 2017, Tim Davis, The University of Aberdeen
    %Creating midpoints and creating element orientations from start and
    %end points of each element
	x = (XBEG + XEND)/2;				% element midpoints
	y = (YBEG + YEND)/2;				% element midpoints
	xe = x;
	ye = y;
	xd = XEND-XBEG;
	yd = YEND-YBEG;
	HalfLength = sqrt(xd.*xd +yd.*yd)/2;			% half-length of each element
 %  In the lines below, B = beta = orientation of element relative to global x-axis
	Beta = atan2(yd,xd);                        % This is the Beta Angle used in the Crouch and Starfield functions
	CosB = xd./(2*HalfLength);					% cosine beta  
    Points=[XBEG,XEND,YBEG,YEND];
    clear XD YD XBEG XEND YBEG YEND
    
    %Now calculating the Normal angle of each element in radians, Presumes
    %normal is 90 deg counter clockwise from orientation segment start-end
    %NormAng is now the outward normal measured from Xaxis counter
    %clockwise to the normal
    NormAng=atan(yd./xd); %
    for i=1:numel(x)
        if xd(i,:)>0
        NormAng(i,:)=NormAng(i,:)+(pi);
        end
    end
    NormAng=NormAng+(pi/2);%adding 90deg
    

    %Check to see if there are duplicate points on the line and calls error
    %if so. By duplicate points the start and end points of a segment are
    %in exactly the same location
    for i=1:NUM
    b111=xd==0;
    b222=yd==0;
    b333=2==(b111+b222);
        if any(b333>0)
        k = find(b333);
        disp 'Locations of duplicate verticies', disp (k);
        error('Duplicate points exist on your line, linear equations will fail, please recreate the line')
        end
    end
end

