function [ Vvaragin,Data ] = AdditionalArgsInVaragin( string,Vvaragin,DefaultVal )
%AdditionalArgsInVaragin: Finding string in the input args varargin (inside 
%               another function), if so we take the arg after use it. We
%               then remove these inputs from the varagin cell. This can be
%               used where a function uses various arguments for plotting
%               etc, additional arguments can also be supplied in these
%               varagin so this function takes these out before the
%               plotting begins.
%
%               Dont use 'nargin' in the function containing this after
%               calling this. Use 'numel(varargin)' as data found by this
%               function is removed from varargin in the outputs. 
% 
%
% usage #1: Gets a flag value defined after 'Flag' in the inputs to a function
% [ varargin,Flag ] = AdditionalArgsInVaragin( 'Flag',varargin )
%               
%
% Arguments: (input)
% Vvaragin      - The varargin arguments inside a function.
%
% string        - A string, e.g. X ='sampling';.
%
% DefaultVal    - If the string is not found the user hasnt input this, we
% 			   therefore jump back onto a default value. 
%
% Arguments: (output)
% Vvaragin      - The input varagin but if data was found inside this has been
%                removed.
%
% Data          - If the argument was found then this is its value, if not we
%                use the default defined in the input arg 'DefaultVal'.
%
%
% Example usage:
%
% % Here 4 arguments are checked for within the varargin introduced into a
% % function. If they are not found these are set to the default value (the
% % last argument). If they are found these are set to the value of the
% % input.
%
% % First we call a function such as (this is its inputs args within the func):
% DrawStressEllipsoidsPrincipal(StressTensor,StrainTensor,X,Y,Z,varargin)
%
% % In this example we call the function as:
% DrawStressEllipsoidsPrincipal(StressTensor,StrainTensor,X,Y,Z,'Scale',4)
%
% %Inside the function is the code: 
%
% [ varargin,Scale ]     = AdditionalArgsInVaragin( 'Scale',    varargin, 1 );
% [ varargin,LineWidth ] = AdditionalArgsInVaragin( 'Sample',   varargin, 2 );
% [ varargin,Points ]    = AdditionalArgsInVaragin( 'Points',   varargin, [] );
% [ ~,Triangles ]        = AdditionalArgsInVaragin( 'Triangles',varargin, [] );
%
% % In this case the argument 'Scale' will be 4 as set in the inputs. 
% % 'Sample' will be its default value of 2 and Points and triangles will
% % be empty variables.
% % The function then uses the variables as it sees fit. 
%
%  Author: Tim Davis
%  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


%Getting the index of the string we are looking for (internal func)
[ Data_indx ] = StringInCellIndex( string,Vvaragin );

%If its empty we exit the function
if isempty(Data_indx)
    Data=DefaultVal;
    %Vvaragin stays the same
else    
	%Getting the data defined after this string	
	Data=Vvaragin{Data_indx+1};
	%Now removing this from the varargin
	Vvaragin(Data_indx:Data_indx+1) = [];
end

function [ Data_indx ] = StringInCellIndex( string,cell )
%StringInCellIndex: Returns the index of a string inside a cell array. 
%
% string - A string, e.g. X ='peter';
%
% cell   - A cell array e.g. G={1,2};

cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));
Data_indx = find(cellfun(cellfind(string),cell));

end

end
