function [varargout] = CreateBlankVars
% CreateBlankVars: Creates empty variables, just call with outputs
%                   variables with titles you would like. Useful if you are
%                   inputting such vars into a function where the input
%                   args are not needed (could be blank arrays) but you
%                   prefer to see variables with names rather than
%                   [],[],[], etc... in the inputs to the function call. 
%               
% usage #1:
% [varargout] = CreateBlankVars
%
% Arguments: (input)
% N/A
%
% Arguments: (output)
% varargout         - A list of output args, the outputs will be vars that
%                    are the same if produced using "A=[];"
%
% Example usage 1:
% %Create 3 variables called A B and C:
% [A,B,C] = CreateBlankVars
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

varargout= cell(1,nargout);

end

