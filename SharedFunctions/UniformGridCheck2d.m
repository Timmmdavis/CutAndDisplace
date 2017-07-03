function [ Uniform ] = UniformGridCheck2d( X,Y )
%UniformGridCheck2d Checks X and Y grid points to see if these have a uniform linear
%spacing. 

%Can be used to check before plotting with different functions or
%for using in interpolation. Returns 1 if true. Can use in if statements 
%i.e.. if uniform spaced points plot with contourf else use scatter. 

%   Copyright 2017, Tim Davis, The University of Aberdeen


if iscolumn(X)
    flagx=1;
    flagy=1;
else    
    %checking that grids are uniformly spaced
    flagx=zeros(1,numel(X(:,1))-1);
    for i=1:numel(X(:,1))-1
        FirstSpace=round(X(1,1)-X(1,2),8);
        spacing=round(X(:,i)-X(:,i+1),8);
        LogEq=spacing==FirstSpace;
        if any(~LogEq)
            flagx(i)=1;
        else
            flagx(i)=0;
        end
    end    
    %Grid data check Y
    flagy=zeros(1,numel(Y(:,1))-1);
    for i=1:numel(Y(:,1))-1
        FirstSpace=round(Y(1,1)-Y(1,2),8);
        spacing=round(Y(:,i)-Y(:,i+1),8);
        LogEq=spacing==FirstSpace;
        if any(~LogEq)
            flagy(i)=1;
        else
            flagy(i)=0;
        end
    end    
end

Uniform= any(~flagx) || any(~flagy); %if the data is a uniform grid


end

