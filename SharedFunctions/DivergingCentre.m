function DivergingCentre( Var )
% DivergingCentre: For a figure colourmapped with a Diverging Cmap use this
%                   to place the 0 value at 0 and not have this some where
%                   at the mid value of the data. 
%               
% usage #1:
% DivergingCentre( Var )
%
% Arguments: (input)
% Var               - This variable is the one used to colour your
%                     colourmap. (Defines the colour values) 
%
% Arguments: (output)
% N/A
%
% Example usage 1:
%
% %Creating a funky div colourmap
% ColourMap1 = [zeros(1, 132), linspace(0, 1, 124)];
% ColourMap2 = [linspace(1, 0, 124), zeros(1, 132)];
% ColourMap = [ColourMap1; ColourMap2; zeros(1, 256)]';
% ColourMap=1-ColourMap;
% %Grabbing from matlabs default data
% Z=peaks;
% %Making sure the data is skewed in the positive axis. 
% Z(Z>0)=Z(Z>0)*2; 
% f1=figure;
% title('Colours centered around data limits')
% surf(Z)
% colormap(ColourMap); colorbar;
% %Copying to a new fig
% f2=copyobj(gcf,0);
% DivergingCentre( Z )
% title('Colours centered around Z=0')
% %Compare the figures, 1st doesnt have the white space at Z~=0. 
% 
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%DivergingCentre For a figure colourmapped with a diverging Cmap use this
%to place the 0 value at 0 and not the mid value

%Turn on the colourbar
colorbar;
%Get the limits of the data, find the maximum.
c = max(abs([min(Var(:)),max(Var(:))]));

if c==0 || isnan(c)
    return
end
%Set the caxis limits to the maximum
clim([-c c]);

end

