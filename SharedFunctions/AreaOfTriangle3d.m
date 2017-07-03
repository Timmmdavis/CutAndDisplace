function [ Area ] = AreaOfTriangle3d( x1,y1,z1,x2,y2,z2,x3,y3,z3 )
%Calculate an area of a tri

%x1,y1,z1 =  xyz location of one corner of the triangle
%x2,y2,z2 =  xyz location of 2nd corner of the triangle
%x3,y3,z3 =  xyz location of 3rd corner of the triangle

%http://math.stackexchange.com/questions/128991/how-to-calculate-area-of-3d-triangle

%To find for multiple tris you call as:
%in my code you would call as:
% Points2=Points(:,2:4);
% for i = 1:numel(Triangles(:,1))
%     %Grabbing the current triangle points for the loop
%     CurrentT1=(Triangles(i,1));
%     CurrentT2=(Triangles(i,2));
%     CurrentT3=(Triangles(i,3));
%     %Getting XYZ list for each vertex on the triangle
%     Pa=Points2(CurrentT1,:);
%     Pb=Points2(CurrentT2,:);
%     Pc=Points2(CurrentT3,:);
%     %Calling func to find areas  
% 	  [ Area ] = AreaOfTriangle3d( Pa(1),Pa(2),Pa(3),Pb(1),Pb(2),Pb(3),Pc(1),Pc(2),Pc(3) );
%     %Appending  
% 	  Data(i)=Area; 
% end	

%   Copyright 2017, Tim Davis, The University of Aberdeen

    a=distance3d(x1,y1,z1,x2,y2,z2);  
    b=distance3d(x2,y2,z2,x3,y3,z3); 
    c=distance3d(x3,y3,z3,x1,y1,z1) ;
    Area = heron(a,b,c)  ;

    function [Area]=heron(a,b,c)  
    s=(a+b+c)/2;   
    Area=(s*(s-a)*(s-b)*(s-c))^0.5;             
    end

    function [d]=distance3d(x1,y1,z1,x2,y2,z2)    
    aa=(x1-x2)^2+(y1-y2)^2+(z1-z2)^2;
    d=aa^0.5;  
    end

end

