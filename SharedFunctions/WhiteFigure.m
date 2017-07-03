function [] = WhiteFigure
%WhiteFigure Makes the current figure turn white (both background and
%outline)
%   Grabs the current figure and makes sure both the background and
%   surrounding box are white. Good for saving the figure. 

%%Example:
%figure;        %create fig
%WhiteFigure;   %make sure its white

%   Copyright 2017, Tim Davis, The University of Aberdeen




set(gcf,'color',[1 1 1])%gcf=Get Current Figure, to set border/figure color
set(gca,'color',[1 1 1])%gca=Get Current Axes, to set axes color

end

