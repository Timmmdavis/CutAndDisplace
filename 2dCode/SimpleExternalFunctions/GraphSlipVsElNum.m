function GraphSlipVsElNum( Dn,Ds )
% GraphSlipVsElNum: Very simple function that draws a plot of slip vs the
%                   element number output from the BEM code.
%                   This plots slips as stair cases.
%
%               
% usage #1:
% GraphSlipVsElNum( Dn,Ds )
%
% Arguments: (input)
%     Dn,Ds        - Two vectors representing the normal (Dn) and shear
%                   (Ds) displacement on each element.
%
%
% Example usage:
%
%  GraphSlipVsElNum( Dn,Ds );
%
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen



NUM=numel(Dn);

%   Plot displacement discontinuity magnitude for all elements
figure,
hold on;
stairs(0:NUM+1,[0;Dn;0],'b'); 
stairs(0:NUM+1,[0;Ds;0],'r');
title('Displacement discontinuity');
xlabel('Element number');
ylabel('Dn (blue) and Ds (red) (m)');
hold off

%Make figure white
WhiteFigure;

end

