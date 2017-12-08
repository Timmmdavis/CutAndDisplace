function DrawS1S2Directions( X,Y,S1dir,S2dir,varargin )
% DrawContourFPlots2d: If S1 and S2 principal directions have been
%                   calculated this function draw these as a quiver plot
%                   Red and Blue.
%               
% usage #1: Just the directions
% DrawS1S2Directions( X,Y,S1dir,S2dir )
%
% usage #2: Fractures also
% DrawS1S2Directions( X,Y,S1dir,S2dir,'P1',P1,'P2',P2,'Scale',0.2 )
%
%
% Arguments: (input)
%     X,Y          - X and Y locations of the data defined in the other
%                   input arguments. Must be a vector.
%
%     S1dir,S2dir  - The direction cosines of the most and least
%                   compressive stress. Locations corresponding to data in
%                   XY. First col should be angle CosAx and 2nd angle
%                   CosAy.
%
% Additional arguments: (input)
%
% Scale,LineWidth   - Properties changing the length and the width of the
%                    resultant lines. Default to '1' and '2' if not set.
%                    MATLAB default width is 0.5.
%
%
%     P1,P2        - Additionally draws the fractures if you want these.
%                   P1 and P2 are the start (P1) and end (P2) points of
%                   each separate element in [x,y] coordiantes. First col
%                   is x.
%
% Example usage:
%  %Assuming we have stress tensors on grid XY:
%  [S1,S2,S1dir,S2dir]=EigCalc2d(Sxx,Syy,Sxy);
%  DrawS1S2Directions( X,Y,S1dir,S2dir )
%
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%Checking if additional arguments have been input into the function:
[ varargin,Scale ]     = AdditionalArgsInVaragin( 'Scale',    varargin,1 );
[ varargin,LineWidth ] = AdditionalArgsInVaragin( 'LineWidth',varargin,2 );
[ varargin,P1 ] = AdditionalArgsInVaragin( 'P1',varargin,[]);
[ ~,P2 ]        = AdditionalArgsInVaragin( 'P2',varargin,[]);

%Extracting parts for visability
S1dirx=S1dir(:,1);
S1diry=S1dir(:,2);
S2dirx=S2dir(:,1);
S2diry=S2dir(:,2);

figure,
hold on
q1  = quiver(X(:),Y(:), S1dirx(:), S1diry(:),LineWidth,'.'); 
q15 = quiver(X(:),Y(:),-S1dirx(:),-S1diry(:),LineWidth,'.');
q2  = quiver(X(:),Y(:), S2dirx(:), S2diry(:),LineWidth,'.');
q25 = quiver(X(:),Y(:),-S2dirx(:),-S2diry(:),LineWidth,'.');
q1.Color = 'red';q15.Color = 'red';q2.Color = 'b';q25.Color = 'b';
q1.AutoScaleFactor=Scale;q15.AutoScaleFactor=Scale;q2.AutoScaleFactor=Scale;q25.AutoScaleFactor=Scale;

chk=exist('P1','var');
if chk==1 
    PlotFracture( P1,P2,'r' )
end

axis equal,
xlabel('x'); ylabel('y'); axis('equal'); title('S1 (red) S2 (blue) directions'); 
hold off

%Make figure white
WhiteFigure

end

