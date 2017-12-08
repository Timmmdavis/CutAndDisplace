function [ C ] = RGB2Colour( v1,v2,v3 )
% RGB2Colour: converts RGB values to colors for MATLAB
%   				 Get RGB values from somewhere like:
%					 http://www.rapidtables.com/web/color/color-wheel.htm
%					 Allowing you to use this to colour lines surfaces etc. 
%               
% usage #1:
% [ C ] = RGB2Colour( v1,v2,v3 )
%
% usage #2:
% RGB=[v1,v2,v3];
% [ C ] = RGB2Colour( RGB )
%
% Arguments: (input)
% v1,v2,v3          - RGB colour values, if you only input v1 this can be a
%                     3*1 matrix i.e. v1=[v1,v2,v3];
%
% Arguments: (output)
% C         		- a 1*3 vector of colour values to be called when plotting.
% 
% Example usage (1):
%
% x = 0:pi/100:2*pi;
% y = sin(x);
% plot(x,y); hold on
% Purp=[191,66,244]; % (this is purple)
% [ C ] = RGB2Colour(Purp(1),Purp(2),Purp(3) ); 
% plot(x+pi,y,'color', C);
% 
% %Example usage (2):
% Green=[0,255,196]; % (this is green)
% [ C ] = RGB2Colour(Green); 
% plot(x+pi*2,y,'color', C);
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

%We assume the user has input this as a 1*3 matrix. 
if nargin==1
	v3=v1(3);
	v2=v1(2);
	v1=v1(1);
end

%Do calc:
C=([v1;v2;v3])/255;

end

