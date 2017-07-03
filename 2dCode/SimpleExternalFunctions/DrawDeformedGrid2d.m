function DrawDeformedGrid2d( X,Y,Ux,Uy,Scl,cmap,ColourVal )
%DrawDeformedGrid2d Draws a wireframe or deformed mesh to give a better
%impression of displacement than quiver plots. Relies on Gridded 2d data. 

%   Copyright 2017, Tim Davis, The University of Aberdeen
sz=size(X); %grabbing the grid size

Xux=X(:)+(Ux(:).*Scl); %adding data to X and Y to deform it
Yuy=Y(:)+(Uy(:).*Scl);

TD=abs(Ux(:))+abs(Uy(:)); %total displacement, we colour for this


TD=reshape(TD,sz); %reshaping back to grid data
Xux=reshape(Xux,sz);
Yuy=reshape(Yuy,sz);

if nargin==6
figure;surf(Xux,Yuy,zeros(size(Yuy)),TD);colormap(cmap);divergingCentre( TD )
else 
ColourVal=reshape(ColourVal,sz);
figure;surf(Xux,Yuy,zeros(size(Yuy)),ColourVal);colormap(cmap) %draw with the imported value
end
az = 0;el = 90;view(az, el); %looking from above
xlabel('x'); ylabel('y'); axis('equal'); title('Displacement')
%camproj('perspective')


end

