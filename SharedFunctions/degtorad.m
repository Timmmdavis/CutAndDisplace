function [ radians ] = degtorad( degrees )
% degtorad: Converts the input argument in degrees to radians, like the
%                   MATLAB function but do not need the toolbox. Note there
%                   is also the MATLAB func: deg2rad() but to save getting
%                   confused between using 'to' and '2' when calling this
%                   it's better to have this function around.
%               
% usage #1:
% [ radians ] = degtorad( degrees )
%
% Arguments: (input)
% degrees           - A vector or matrix with values in degrees. 
%
% Arguments: (output)
% radians           - A vector or matrix with the input arguments converted
%                    to radians.
%
% Example usage 1:
% [ SixtySixInRad ] = degtorad( 66 )
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%Function:
radians=degrees.*(pi./180);

end

