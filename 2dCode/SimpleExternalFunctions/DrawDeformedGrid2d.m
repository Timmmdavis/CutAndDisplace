function DrawDeformedGrid2d( X,Y,Ux,Uy,cmap,ColourVal,varargin )
% DrawDeformedGrid2d: Draws a wireframe or deformed mesh to give a better
%                   impression of displacement than quiver plots. Relies
%                   on Gridded 2d data.
%                   
%                   Used for 2D plots in: Davis, T., Healy, D., Bubeck, A.
%                   and Walker, R., 2017. Stress concentrations around
%                   voids in three dimensions: The roots of failure.
%                   Journal of Structural Geology, 102, pp.193-207.
%               
% usage #1:
% DrawDeformedGrid2d( X,Y,Ux,Uy,Scale,cmap,ColourVal )
%
% usage #2: No colourmap defined, defaults to grey. 
% DrawDeformedGrid2d( X,Y,Ux,Uy,Scale,[],ColourVal )
%
% Arguments: (input)
%     X,Y          - X and Y locations of the Ux and Uy data defined in the
%                   other input arguments. Must be on a grid.
%
%     Ux,Uy        - Cartesian displacement components at points XY. Must
%                   match the size and shape of XY.
%
%     cmap         - A colourmap that MATLAB can use. See func
%                   "colormap_cpt.m" to produce one. 
%
%     ColourVal    - Colours the deformation by this value. Must be the
%                   size and shape of the other inputs.
%
% Additional arguments: (input)
%     Scale        - Scales the deformation for viewing. If the plot
%                   overlaps itself reduce the scale. If you see no
%                   deformation increase this.
%
% Example usage:
%
%  [X,Y]=meshgrid(-2:0.1:2,-2:0.1:2);
%  Ux=Y/5;Uy=X/5;
%  maxdisp=sqrt(Ux.^2+Uy.^2);
%  cmap=flipud(colormap(gray));
%  DrawDeformedGrid2d( X,Y,Ux,Uy,cmap,maxdisp,'Scale',10 );
%
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

[ ~,Scale ]  = AdditionalArgsInVaragin( 'Scale', varargin,1 );

sz=size(X); %grabbing the grid size

Xux=X(:)+(Ux(:).*Scale); %adding data to X and Y to deform it
Yuy=Y(:)+(Uy(:).*Scale);

TD=abs(Ux(:))+abs(Uy(:)); %total displacement, we colour for this
TD=reshape(TD,sz); %reshaping back to grid data
Xux=reshape(Xux,sz);
Yuy=reshape(Yuy,sz);

figure;
if nargin==5 %No Data to colour for supplied
    surf(Xux,Yuy,zeros(size(Yuy)),TD);
    colormap(cmap);
    DivergingCentre( TD )
else 
    surf(Xux,Yuy,zeros(size(Yuy)),ColourVal);
    if isempty(cmap) 
        colormap(gray); 
    else 
        colormap(cmap); 
    end %draw with the imported value
end
az = 0;el = 90;view(az, el); %looking from above
xlabel('x'); ylabel('y'); axis('equal'); title('Displacement')
%camproj('perspective')

%Make figure white
WhiteFigure


end

