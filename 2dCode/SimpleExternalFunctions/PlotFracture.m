function PlotFracture( P1,P2,colour )
%PlotFracture Plots the loaded 2d fractures and the normals of each
%element. 

hold on    
line([P1(:,1)';P2(:,1)'],[P1(:,2)';P2(:,2)'],'color',colour)

%Make figure white
WhiteFigure;

end

