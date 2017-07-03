function divergingCentre( Var )
%divergingCENTRE For a figure colourmapped with a diverging Cmap use this
%to place the 0 value at 0 and not the mid value

%Usage:
% figure;                                                   %create blank fig  
% scatter(X(:),Y(:),15,Sxy(:))                              %put values in fig
% colormap(cmap2)                                           %load a diverging map as the current Cmap
% xlabel('x'); ylabel('y'); axis('equal'); title('Sxy');    %put labels etc on fig
% divergingCentre( Sxy )                                    %Use this func to control the cmap

%   Copyright 2017, Tim Davis, The University of Aberdeen

colorbar;
c = max(abs([min(Var(:)),max(Var(:))]));
caxis([-c c]);



end

