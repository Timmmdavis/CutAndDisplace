
%Clear all
c

e=linspace(-1,1,50);
 N1=((-1/3)*e+(2/3)*e.^2);
 N2=(1-(4/3)*e.^2);
 N3=((1/3)*e+(2/3)*e.^2);

%Locations on element.
Left=-(sqrt(3))/2;
Mid=0;
Right=(sqrt(3))/2;

%Different shape functions of disp. 
N1p=((-1/3)*Mid+(2/3)*Mid.^2);
N2p=(1-(4/3)*Mid.^2);
N3p=((1/3)*Mid+(2/3)*Mid.^2);

scatter(e,N1);
hold on
scatter(Mid,N1p);