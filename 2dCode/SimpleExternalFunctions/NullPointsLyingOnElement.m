function [ X,Y  ] = NullPointsLyingOnElement( X,Y,Points,Midpoint )
%NullPointsLyingOnElement Finds any observation points lying on the element
%and puts these to NAN
%WHAT ABOUT MIDPOINTS?

XBEG=Points(:,1);
XEND=Points(:,2);
YBEG=Points(:,3);
YEND=Points(:,4);

NUM=numel(XBEG);

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

