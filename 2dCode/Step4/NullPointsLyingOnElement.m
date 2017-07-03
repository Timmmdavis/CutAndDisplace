%   Copyright 2017, Tim Davis, The University of Aberdeen

function [ X,Y  ] = NullPointsLyingOnElement( X,Y,XBEG,YBEG,XEND,YEND,NUM )
%NullPointsLyingOnElement Finds any observation points lying on the element
%and puts these to NAN
for i=1:NUM
    b111=X==XBEG(i);
    b222=Y==YBEG(i);
    b333=2==(b111+b222);
    Y(b333)=nan; X(b333)=nan;
    
    b111=X==XEND(i);
    b222=Y==YEND(i);
    b333=2==(b111+b222);
    Y(b333)=nan; X(b333)=nan;
end    

end

