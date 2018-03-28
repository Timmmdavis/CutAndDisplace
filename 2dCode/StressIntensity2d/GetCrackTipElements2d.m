function [EndElsLoc] = GetCrackTipElements2d(P1,P2)
% GetCrackTipElements2d: Gets the end points of 2D elements that are 'Free'. 
%                   I.e. not connected to another element. Essentially just
%                   want to find rows in P1 or P2 that are unique. The
%                   output is a flag where values above 0 are end elements.
%                   If the value is 1 then the end point is in P1 and if 2
%                   in P2.
%                   Created to find elements where stress intensity is
%                   computed. 
%               
% usage:
% [EndElsLoc] = GetCrackTipElements2d(P1,P2)
%
% Arguments: (input)
% P1,P2              - The end points of each line segment. [X,Y] 
%                     Arranged so the row index's correspond to HalfLength
%                     etc. 
%
% Arguments: (output)
% EndElsLoc          - A vector the same length as P1. Acts as a flag to
%                      say if an element is an endpoint. The number in the
%                      flag (1 or 2) says if its P1 or P2 that is the
%                      actual end/free point. 
%
% Example usage:
%
% x=linspace(-0.5,0.5,15);    
% y=zeros(1,numel(x)); 
% Pointsxy=[x;y]';
% mystruct.line1=(1:(length(Pointsxy(:,1))));
% [MidPoint,HalfLength,P1,P2,LineNormalVector,Lines,Points,Fdisp  ]...
%  = CreateElements2d( Pointsxy,mystruct,[] );
% [EndElsLoc] = GetEndElements2d(P1,P2);
% Flg=EndElsLoc>0;
% PlotFracture( P1,P2,'r' ); hold on
% PlotFracture( P1(Flg,:),P2(Flg,:),'b' );
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University



%Check if one is inside the other. 0's here are end elements. 
EndElsP1 = ismember(P1,P2,'rows');
EndElsP2 = ismember(P2,P1,'rows');
%Create flag
EndElsLoc=zeros(size(EndElsP1));
%Add end elements that are in P1 (set as 1's)
EndElsLoc=EndElsLoc+(EndElsP1==0); 
%Add end elements that are in P2 (set as 2's)
EndElsLoc=EndElsLoc+(EndElsP2==0)*2; 

%Just checking elements with normals pointing the other direction do not
%mess this up. 
if numel(unique(P1,'rows'))~=numel(P1)
    error('P1 has two elements with the same end point')
end
if numel(unique(P2,'rows'))~=numel(P2)
    error('P2 has two elements with the same end point')
end

end

