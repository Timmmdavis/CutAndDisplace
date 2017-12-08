function [ degrees ] = radtodeg( radians )
% radtodeg: Converts radians to degrees, like the
%                   MATLAB function but do not need the toolbox. Note there
%                   is also the MATLAB func: rad2deg() but to save getting
%                   confused between using 'to' and '2' when calling this
%                   it's better to have this function around.
%
% usage #1: 
% [ degrees ] = radtodeg( radians )
%                 
% Arguments: (input)
% radians       - The data in radians that will be converted to degrees.
%
% Arguments: (input)
% degrees       - The converted data in degrees (shape maintained). 
%
% Example usage:
%
% radians = [1 2; 3 2.2];
% [ degrees ] = radtodeg( radians )
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%Function:
degrees=radians.*(180./pi);

end

