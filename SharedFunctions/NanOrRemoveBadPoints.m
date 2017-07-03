function [varargout] = NanOrRemoveBadPoints( Value,type,absflag,varargin )
%NanOrRemoveBadPoints 
%Value filter value above below or = to this is set as bad
%Type set if you want above below or = to this and nans or remvoving this
%data. 
%absflag, if 1 use abs values. 
%varargin the number of arguments you are passing in. 

%   Copyright 2017, Tim Davis, The University of Aberdeen

%%%%%Removing will turn data to col vec
%type=1 values above Value removed
%type=2 values below Value removed
%type=3 values same as Value removed
%type=4 values flagged by logical 'Value' changed removed;

%%%%%Nan will keep values the same. Octave contourf has some trouble drawing nans
%type=5 values above Value changed to Nan
%type=6 values below Value changed to Nan
%type=7values same as Value changed to Nan
%type=8 values flagged by logical 'Value' changed to Nan;
    
for i=1:nargin-3
    
    
    Orig=varargin{i};
    if absflag==0
        CurVal=Orig;
    elseif absflag==1
        CurVal=abs(Orig);
    end
        
    if (type ==1 ||  type ==5) %1 or 5
    bad=CurVal>Value;
    elseif (type ==2 ||  type ==6) 
    bad=CurVal<Value;    
    elseif (type ==3 ||  type ==7)   
    bad=CurVal<Value;    
    elseif (type ==4 ||  type ==8)   
    bad=Value;    
    end
    
    if type>4 %change to nan
    Orig(bad)=nan;
    else      %remove data  
    Orig(bad)=[];        
    end
    varargin{i}=Orig;
end


varargout=varargin;


end

