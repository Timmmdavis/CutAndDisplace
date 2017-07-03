%Annulus
function [rOvb,SrrNorm,SttNorm]=AnnulusEq_Func(G1,v1,G2,v2);

%   Copyright 2017, Tim Davis, The University of Aberdeen

%%%%
%Parameters
%%%%
a=1;        %Hole radius
b=2;        %Annulus outer edge to centre radius
Rad_a_r_b=linspace(a,b,50);  %   Vector of observation points lying r from the centre in the annulus
Rad_b_r=linspace(b,3,50);    %   Vector of observation points lying r from the centre in the matrix

% G1=1;       %Elastic properties of Annulus 
% v1=0.25; 
% 
% G2=0.5;     %Elastic properties of Matrix   
% v2=0.25;

p=-10^-3;    %Internal hole pressure  %neg as its compressive  




%%%%
%Equations for PPrime (assigned pr in script)
%%%%
one=2*(1-v1)*p*(a^2)/(b^2);
two=2*(1-v1)+(G1/G2-1)*(1-(a^2)/(b^2));
pr= one/two; %C&S 7.5.24

%%%%
%Equations for calculating Srr (Eq1 crouch and star 7.5.23)
%%%%
one=(1/(1-(a^2)/(b^2)));
two=(((p*(a^2))/(b^2)-pr)-(p-pr)*(a^2)./(Rad_a_r_b.^2));
Srr_a_r_b=one*two;
Srr_b_r=-pr*(b^2)./(Rad_b_r.^2);

%%%%
%Equations for calculating Stt (Eq1 crouch and star 7.5.23)
%%%%
one=(1/(1-(a^2)/(b^2)));
two=(((p*(a^2))/(b^2)-pr)+(p-pr)*(a^2)./(Rad_a_r_b.^2));
Stt_a_r_b=one*two;
Stt_b_r=+pr*(b^2)./(Rad_b_r.^2);

%%%%
%Normalizing
%%%%
SrrNorm_a_r_b=Srr_a_r_b/p;
SrrNorm_b_r=Srr_b_r/p;

SttNorm_a_r_b=Stt_a_r_b/p;
SttNorm_b_r=Stt_b_r/p;

rOvb_a_r_b=Rad_a_r_b/b;
rOvb_b_r=Rad_b_r/b;

%%%%
%Appending
%%%%
SrrNorm=[SrrNorm_a_r_b,SrrNorm_b_r];
SttNorm=[SttNorm_a_r_b,SttNorm_b_r];
rOvb=[rOvb_a_r_b,rOvb_b_r];

%%%%
%Plotting
%%%%
% plot(rOvb,SrrNorm);
% hold on
% plot(rOvb,SttNorm);
% title('Srr (blue) Stt (orange)'), xlabel('rad/b')
% ylabel('SrrSttNorm');
