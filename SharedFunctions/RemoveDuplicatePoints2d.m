function [ X,Y ] = RemoveDuplicatePoints2d( X,Y )
%RemoveDuplicatePoints2d: Finds and removes any duplicate (coincident)
%               points for 2D X Y data. 
% usage #1:
% [ X,Y ] = RemoveDuplicatePoints2d( X,Y )
%
% Arguments: (input)
% X             - list of x points (vect)
%
% Y             - list of y points (vect)
%
% Arguments: (output)
% X             - list of x points (vect) with duplicate points removed
%
% Y             - list of y points (vect) with duplicate points removed
% 
%
% Example usage:
%
% %Create a square of points (first and last duplicate)
% x = [0 0 1 1 0]; 
% y = [0 1 1 0 0]; 
% %Removes the last point that is a duplicate
% [ x,y ] = RemoveDuplicatePoints2d( x,y )
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%Repeat data
Xlst=X;
Ylst=Y;
%Start loop searching for similar data. 
for i=1:numel(X)
    DupPnts=(Xlst(i)==X)+(Ylst(i)==Y)==2;
    %If there are duplcate points
    if sum(DupPnts)>1
        %Get indx of these
        indx = find(DupPnts); 
        %duplicate point (nan the latter points)
        X(indx(2:end))=nan; Y(indx(2:end))=nan;
    end
end    
%Now removing these
X(isnan(X))=[];
Y(isnan(Y))=[];


end

