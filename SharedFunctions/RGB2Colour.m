function [ C ] = RGB2Colour( v1,v2,v3 )
%RGB2Colour converts RGB values to colors for MATLAB
%   Get RGB values from somewhere like: http://www.rapidtables.com/web/color/color-wheel.htm
%   Values like :191,66,244 (this is purple)
%   put these into this func: [ C ] = RGB2Colour( 191,66,244 ); 
%   Then when plotting call: plot(x,y,'color', C);

%   Copyright 2017, Tim Davis, The University of Aberdeen

%C=normc([v1;v2;v3]);
C=([v1;v2;v3])/255;

end

