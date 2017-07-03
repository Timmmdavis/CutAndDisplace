function GraphSlipVsElNum( NUM,TensileDisp,ShearDisp )
%GraphSlipVsElNum draws a plot of slip vs the element number. 
%Plot numerical slips as stair cases

%   Copyright 2017, Tim Davis, The University of Aberdeen


%   Plot displacement discontinuity magnitude for all elements
figure, stairs(0:NUM+1,[0;TensileDisp;0],'b'); hold on; stairs(0:NUM+1,[0;ShearDisp;0],'r');
title('displacement discontinuity'), xlabel('element number')
ylabel('Dn (blue) and Ds (red) (m)'); hold off

end

