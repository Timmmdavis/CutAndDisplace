function  ChangeFontSizes( axistickssz,titlesz )
%Change Font sizes in current figure Summary of this function goes here
%  titlesz - size of the fig titles
%  titlesz - size of the axis ticks 
%just bring these in as numbers
%call like: ChangeFontSizes( 15,20 )

%   Copyright 2017, Tim Davis, The University of Aberdeen

set(findall(gcf,'type','axes'),'fontsize',axistickssz)
set(findall(gcf,'type','text'),'fontSize',titlesz) 

set(gcf,'color','w'); %sets background white

%changing colorbar font only if it exists
figure_handle=gcf;
cbar_handle = findobj(figure_handle,'tag','Colorbar');
if ~sum(size(cbar_handle))==0
handle=colorbar;
set(handle,'fontsize',axistickssz);
end


end

