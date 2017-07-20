function [StrikeSlipDisp,DipSlipDisp,TensileSlipDisp,Pxx,Pyy,Pzz,...
Pxy,Pxz,Pyz]=SlipCalculator3d(MidPoint,Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,...
Tnn_i,Tss_i,Tds_i,mu,lambda,nu,P1,P2,P3,halfspace,FaceNormalVector,Fdisp,strain,Mu,Sf,Option,Triangles,Points,cmap)
%SlipCalculator3d Calculates the 3 slip components on the imported fault
%surface from boundary condition defined by the user.
%
%   Influence Coefficient functions from:
%   Nikkhoo, M. and Walter, T.R., 2015. Triangular dislocation: an analytical
%   , artefact-free solution. Geophysical Journal International, 201(2), 
%   pp.1117-1139.
%
%   INPUTS
%   nu is the Poisson's ratio
%   Lambda is the Lamé's parameter
%   mu is the shear modulus
%   halfspace defines if we work out the coefficientsin a half or whole
%   space
%   FaceNormalVector is cosines of each tris normal
%   Pxx,Pyy,Pxy,Tnn,Tss,Tds are the defined stresses or tractions at the
%   elements, either defined for each element or just a single value that
%   will then be added to all elements.
%   strain is a flag that means remote stress are actually strains, these
%   are converted to stresses before this continues. 
%   Fdisp, flag of elements that have fixed displacements. Displacement
%   coefficients are created if this exists and the system of linear
%   equations changes. 
%   Mu, Coefficient of friction
%   Sf, Sliding friction
%   Option - which option the user is using. 
    %B=under remote stress, fault can interpenetrate
    %C=friction solver as in ritz 2012
    %D=remote stress, opening components not calculated. Fault satisfys
    %boundary conditions by shearing only.
    %E=tractions defined on body, Ie pressure in a magma chamber
    %F=Inhomogeneous
    
%   Copyright 2017, Tim Davis, The University of Aberdeen
    
% If Defined slip is as a uniform remote strain, checking options that have
% a uniform 'stress' defined
if Option=='B' || Option=='C' || Option=='D' || Option=='F'

    %for the imported stress making sure these are col vectors
    [ Pxx ] = RowVecToCol( Pxx );
    [ Pyy ] = RowVecToCol( Pyy );
    [ Pzz ] = RowVecToCol( Pzz );
    [ Pxy ] = RowVecToCol( Pxy );
    [ Pxz ] = RowVecToCol( Pxz );
    [ Pyz ] = RowVecToCol( Pyz );

    if strain==1
        %Converting the input strain boundary condition on every faults midpoint to stress
        %Hooke's law, Equation 7.131 and 7.132 in David Pollards Book.  
        [ Pxx,Pyy,Pzz,Pxy,Pxz,Pyz ] = HookesLaw3dStrain2Stress( Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,lambda,mu );
    end  

end

%Function to check the boundary stress tensors related to Z(y in 2d) components are 0 at the
%halfspace surface
if halfspace==1
    HalfSpaceBoundaryConditionsCheck3d(MidPoint(:,3),Pzz,Pyz,Pxz)
else
end
    
    %Checking we can overcome the frictional resistance on the fault. No
    %point waiting for this otherwise
if Option=='C'     
    
    [ Tnn,Tds,Tss ] = CalculateNormalAndShearTractions3d( FaceNormalVector,Pxx,Pyy,Pzz,Pxy,Pxz,Pyz );
    FricRes=Tnn.*Mu; %Negative is compression
    if all(-abs(Tss)>FricRes) && all(-abs(Tds)>FricRes)
        PlotSlipDistribution3d(Triangles,Points,cmap,FricRes,Tss,Tds)
        disp('Your frictional fault will not slip, quit (ctrl+c)')
		disp('TurnedOfPauseForLoopsThoughPUTBACKON!')
        %pause        
    end
    
end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Creating influence matrices using a function. Big square matrices of
    %how much every element slipping by 1 effects traction on every other element
    %If displacement if fixed on any elements we also disp influence matrices 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
if Option~='F'  %any option but inhomo   
    
    [ DnTn,DnTss,DnTds,DssTn,DssTss,DssTds,DdsTn,DdsTss,DdsTds,Dn_dx,Dn_dy,Dn_dz,Dss_dx,Dss_dy,Dss_dz,Dds_dx,Dds_dy,Dds_dz,NUM,...
    StrikeSlipCosine,DipSlipCosine] = CalculateInfluenceMatrices3d(MidPoint,P1,P2,P3,mu,lambda,FaceNormalVector,halfspace,nu,Fdisp);
    
else %Inhomo calc
    
    [DnTnE1FB,DnTssE1FB,DnTdsE1FB,DssTnE1FB,DssTssE1FB,DssTdsE1FB,DdsTnE1FB,DdsTssE1FB,DdsTdsE1FB,...
    DnTnE1IF,DnTssE1IF,DnTdsE1IF,DssTnE1IF,DssTssE1IF,DssTdsE1IF,DdsTnE1IF,DdsTssE1IF,DdsTdsE1IF,...
    Dn_dxE1FB,Dn_dyE1FB,Dn_dzE1FB,Dss_dxE1FB,Dss_dyE1FB,Dss_dzE1FB,Dds_dxE1FB,Dds_dyE1FB,Dds_dzE1FB,...
    Dn_dxE1IF,Dn_dyE1IF,Dn_dzE1IF,Dss_dxE1IF,Dss_dyE1IF,Dss_dzE1IF,Dds_dxE1IF,Dds_dyE1IF,Dds_dzE1IF,...
    DnTnE2FB,DnTssE2FB,DnTdsE2FB,DssTnE2FB,DssTssE2FB,DssTdsE2FB,DdsTnE2FB,DdsTssE2FB,DdsTdsE2FB,...
    DnTnE2IF,DnTssE2IF,DnTdsE2IF,DssTnE2IF,DssTssE2IF,DssTdsE2IF,DdsTnE2IF,DdsTssE2IF,DdsTdsE2IF,...
    Dn_dxE2FB,Dn_dyE2FB,Dn_dzE2FB,Dss_dxE2FB,Dss_dyE2FB,Dss_dzE2FB,Dds_dxE2FB,Dds_dyE2FB,Dds_dzE2FB,...
    Dn_dxE2IF,Dn_dyE2IF,Dn_dzE2IF,Dss_dxE2IF,Dss_dyE2IF,Dss_dzE2IF,Dds_dxE2IF,Dds_dyE2IF,Dds_dzE2IF,...
    NUME1,NUME2,nuE1,muE1,nuE2,muE2,...
    FreeBoundary1,FreeBoundary2,FdispE1,FdispE2,FreeBoundaries,FixedDisps]...
    = InhomogeneousInfluenceBuilder3d...
    (MidPoint,P1,P2,P3,mu,lambda,FaceNormalVector,halfspace,nu,Fdisp);
    
end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Traction influence matrices should now exist. 
    %Now the defined pressures, Dss,Dds and Dn are reshaped. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Option=='B' || Option=='C' || Option=='D' || Option=='F' 
    
    if Option=='F'   

        FaceNormalVector=FaceNormalVector([FreeBoundaries,FreeBoundaries,FreeBoundaries]); 
        FaceNormalVector=reshape(FaceNormalVector,[],3);  %after grabbing with flag turns to Cv, need to reshape  
    end     
    
    [ Tnn,Tds,Tss ] = CalculateNormalAndShearTractions3d( FaceNormalVector,Pxx,Pyy,Pzz,Pxy,Pxz,Pyz );
    
end

    %for option E traction has been already defined
if Option=='E' || ~isempty(Tnn_i)&&~isempty(Tds_i)&&~isempty(Tss_i)
    
    %If Tractions do not exist yet (i.e. option E) then we create blank ones
    %to append the user defined ones too. 
    Tnn=zeros(NUM,1);
    Tss=zeros(NUM,1);
    Tds=zeros(NUM,1);
    
    %creates a column of the stress on every triangles centre point with 
    %the pressure direction defined by the slip type in the name. 
    Tnn     = Tnn+Tnn_i;    
    Tss     = Tss+Tss_i;
    Tds     = Tds+Tds_i;

end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Creating large arrays from the influence matrices and the traction matrices.
    %performing the system of linear equations to find the slip due to the boundary conditions. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Quick loop to force any fixed triangles to stick
if Option~='F'  
    
    if any(Fdisp)==1

    %Now exchanging stress coeff rows for displacement ones
    [ DnTn,DnTss,DnTds,DssTn,DssTss,DssTds,DdsTn,DdsTss,DdsTds,NUM,Tnn,Tss,Tds ]...
    = FixingDisp_InfRowSwap3d( DnTn,DnTss,DnTds,DssTn,DssTss,DssTds,DdsTn,...
    DdsTss,DdsTds,Dn_dx,Dn_dy,Dn_dz,Dss_dx,Dss_dy,Dss_dz,Dds_dx,Dds_dy,Dds_dz,NUM,Tnn,Tss,Tds,Fdisp);
    clear Dn_dx Dn_dy Dn_dz Dss_dx Dss_dy Dss_dz Dds_dx Dds_dy Dds_dz
    end
    
else %We need to do this for both elastics  
    
    if any(FdispE1)==1    
    [DnTnE1FB,DnTssE1FB,DnTdsE1FB,DssTnE1FB,DssTssE1FB,DssTdsE1FB,DdsTnE1FB,DdsTssE1FB,DdsTdsE1FB,NUME1,null,null,null ]...
    = FixingDisp_InfRowSwap3d(DnTnE1FB,DnTssE1FB,DnTdsE1FB,DssTnE1FB,DssTssE1FB,DssTdsE1FB,DdsTnE1FB,DdsTssE1FB,DdsTdsE1FB,...
    Dn_dxE1FB,Dn_dyE1FB,Dn_dzE1FB,Dss_dxE1FB,Dss_dyE1FB,Dss_dzE1FB,Dds_dxE1FB,Dds_dyE1FB,Dds_dzE1FB,NUME1,[],[],[],FdispE1);
    % % %Now removing any fixed disp element cols from the interface elements  
    [DnTnE1IF,DnTssE1IF,DnTdsE1IF,DssTnE1IF,DssTssE1IF,DssTdsE1IF,DdsTnE1IF,DdsTssE1IF,DdsTdsE1IF]...
    = ExtractData(FdispE1,4,DnTnE1IF,DnTssE1IF,DnTdsE1IF,DssTnE1IF,DssTssE1IF,DssTdsE1IF,DdsTnE1IF,DdsTssE1IF,DdsTdsE1IF );
    [Dn_dxE1IF,Dn_dyE1IF,Dn_dzE1IF,Dss_dxE1IF,Dss_dyE1IF,Dss_dzE1IF,Dds_dxE1IF,Dds_dyE1IF,Dds_dzE1IF]...
    = ExtractData(FdispE1,4,Dn_dxE1IF,Dn_dyE1IF,Dn_dzE1IF,Dss_dxE1IF,Dss_dyE1IF,Dss_dzE1IF,Dds_dxE1IF,Dds_dyE1IF,Dds_dzE1IF);
    end
    %E2    
    if any(FdispE2)==1
        
    %Need to Flip So the fixed disps are the top part    
    [DnTnE2FB,DnTssE2FB,DnTdsE2FB,DssTnE2FB,DssTssE2FB,DssTdsE2FB,DdsTnE2FB,DdsTssE2FB,DdsTdsE2FB,FdispE2]...
        =FlipColVecs(DnTnE2FB,DnTssE2FB,DnTdsE2FB,DssTnE2FB,DssTssE2FB,DssTdsE2FB,DdsTnE2FB,DdsTssE2FB,DdsTdsE2FB,FdispE2);
            
    [DnTnE2FB,DnTssE2FB,DnTdsE2FB,DssTnE2FB,DssTssE2FB,DssTdsE2FB,DdsTnE2FB,DdsTssE2FB,DdsTdsE2FB,NUME2,null,null,null ]...
    = FixingDisp_InfRowSwap3d(DnTnE2FB,DnTssE2FB,DnTdsE2FB,DssTnE2FB,DssTssE2FB,DssTdsE2FB,DdsTnE2FB,DdsTssE2FB,DdsTdsE2FB,...
    Dn_dxE2FB,Dn_dyE2FB,Dn_dzE2FB,Dss_dxE2FB,Dss_dyE2FB,Dss_dzE2FB,Dds_dxE2FB,Dds_dyE2FB,Dds_dzE2FB,NUME2,[],[],[],FdispE2);
    
    %Need to Flip back
    [DnTnE2FB,DnTssE2FB,DnTdsE2FB,DssTnE2FB,DssTssE2FB,DssTdsE2FB,DdsTnE2FB,DdsTssE2FB,DdsTdsE2FB,FdispE2]...
        =FlipColVecs(DnTnE2FB,DnTssE2FB,DnTdsE2FB,DssTnE2FB,DssTssE2FB,DssTdsE2FB,DdsTnE2FB,DdsTssE2FB,DdsTdsE2FB,FdispE2);
            
    % % %Now removing any fixed disp element cols from the interface elements  
    [DnTnE2IF,DnTssE2IF,DnTdsE2IF,DssTnE2IF,DssTssE2IF,DssTdsE2IF,DdsTnE2IF,DdsTssE2IF,DdsTdsE2IF]...
    = ExtractData(FdispE2,4,DnTnE2IF,DnTssE2IF,DnTdsE2IF,DssTnE2IF,DssTssE2IF,DssTdsE2IF,DdsTnE2IF,DdsTssE2IF,DdsTdsE2IF );
    [Dn_dxE2IF,Dn_dyE2IF,Dn_dzE2IF,Dss_dxE2IF,Dss_dyE2IF,Dss_dzE2IF,Dds_dxE2IF,Dds_dyE2IF,Dds_dzE2IF]...
    = ExtractData(FdispE2,4,Dn_dxE2IF,Dn_dyE2IF,Dn_dzE2IF,Dss_dxE2IF,Dss_dyE2IF,Dss_dzE2IF,Dds_dxE2IF,Dds_dyE2IF,Dds_dzE2IF);       
    end
    %Need to do this for Tn and Ts 
    FD=logical(FixedDisps);    
    Tnn(FD) = 0;    
    Tss(FD) = 0;
    Tds(FD) = 0;
    
end       
    
if Option=='B' || Option=='E' || Option=='C'
    Ats = [-DnTn,  -DssTn, -DdsTn]; %geological convention from coeff matrix builder, we want engin
    Ass = [-DnTss, -DssTss,-DdsTss];
    Ads = [-DnTds, -DssTds,-DdsTds];
    clear DnTn DssTn DdsTn DnTssDssTss DdsTss DnTds DssTds DdsTds
    A= [Ats;Ass;Ads];  %Concatenate ready for equation 
    clear Ats Ass Ads

elseif Option=='D'                                                
    Ats = [-DssTn, -DdsTn]; %geological convention from coeff matrix builder, we want engin
    Ass = [-DssTss,-DdsTss];
    Ads = [-DssTds,-DdsTds];
    clear DnTn DssTn DdsTn DnTssDssTss DdsTss DnTds DssTds DdsTds
    A= [Ats;Ass;Ads];  %Concatenate ready for equation 
    clear Ats Ass Ads
    
elseif Option=='F'         
    zerinfA=zeros(size(DnTnE1FB,1),size(Dn_dzE2IF,2)); %Sz=1=(Rows down),2=(Cols across)
    zerinfB=zeros(size(DnTnE2FB,1),size(Dn_dzE1IF,2)); 
    
    A=[[DnTnE1IF,  DssTnE1IF,  DdsTnE1IF,   -DnTnE2IF,     -DssTnE2IF,    -DdsTnE2IF  ];
       [DnTssE1IF, DssTssE1IF, DdsTssE1IF,  -DnTssE2IF,    -DssTssE2IF,   -DdsTssE2IF ];
       [DnTdsE1IF, DssTdsE1IF, DdsTdsE1IF,  -DnTdsE2IF,    -DssTdsE2IF,   -DdsTdsE2IF ];
       [Dn_dxE1IF, Dss_dxE1IF,  Dds_dxE1IF,   -Dn_dxE2IF,    -Dss_dxE2IF,    -Dds_dxE2IF  ];
       [Dn_dyE1IF, Dss_dyE1IF,  Dds_dyE1IF,   -Dn_dyE2IF,    -Dss_dyE2IF,    -Dds_dyE2IF  ];
       [Dn_dzE1IF, Dss_dzE1IF,  Dds_dzE1IF,   -Dn_dzE2IF,    -Dss_dzE2IF,    -Dds_dzE2IF  ];
       [DnTnE1FB,  DssTnE1FB,  DdsTnE1FB,   zerinfA,         zerinfA,         zerinfA ];
       [DnTssE1FB, DssTssE1FB, DdsTssE1FB,  zerinfA,         zerinfA,         zerinfA ];
       [DnTdsE1FB, DssTdsE1FB, DdsTdsE1FB,  zerinfA,         zerinfA,         zerinfA ];    
       [zerinfB,    zerinfB,   zerinfB,     DnTnE2FB,       DssTnE2FB,     DdsTnE2FB  ]; %stress infs, only the influence on freeboundary
       [zerinfB,    zerinfB,   zerinfB,     DnTssE2FB,      DssTssE2FB,    DdsTssE2FB ];
       [zerinfB,    zerinfB,   zerinfB,     DnTdsE2FB,      DssTdsE2FB,    DdsTdsE2FB]];
    
    clear DnTnE1IF DnTssE1IF DnTdsE1IF DssTnE1IF DssTssE1IF DssTdsE1IF DdsTnE1IF DdsTssE1IF DdsTdsE1IF zerinfSm
    clear Dn_dxE2IF Dn_dyE2IF Dss_dxE2IF Dss_dyE2IF Dss_dzE2IF Dds_dxE2IF Dds_dyE2IF Dds_dzE2IF
    clear DnTnE1FB DnTssE1FB DnTdsE1FB DssTnE1FB DssTssE1FB DssTdsE1FB DdsTnE1FB DdsTssE1FB DdsTdsE1FB
    clear DnTnE2FB DnTssE2FB DnTdsE2FB DssTnE2FB DssTssE2FB DssTdsE2FB DdsTnE2FB DdsTssE2FB DdsTdsE2FB
   
end

%    [ A ] = InfMatrixCheck( A );  
     disp('inf check off, line 243 in *slipcalculator3d*, put back on'); 
    
if Option=='F'
    zer=zeros(size(Dn_dzE2IF,1),1);     clear Dn_dzE2IF
    contvec=[zer;zer;zer;zer;zer;zer];
    B=-[contvec;...
        Tnn(FreeBoundary1);Tss(FreeBoundary1);Tds(FreeBoundary1);...
        Tnn(FreeBoundary2);Tss(FreeBoundary2);Tds(FreeBoundary2)];
else    
    %Creating vector of tractions at every element
    B= [Tnn;-Tss;Tds];
end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %The linear equations: 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %  Solve the matrix system [A][D]=[B] to get the displacement discontinuity vector D
    %If we have fixed displacements the matrix is no longer square, in this
    %case \ division is faster on sparse matrices 
    if any(Fdisp)==1 && isa(A,'double') %sparse doesn't work for singles
    %D = A\B; disp('sparse off, line 260 SlipCalc3d')
    D = sparse(A)\B;
    %If not we are using a dense square matrix, this is faster if square as we
    %do not need to spend time allocating this as sparse. 
    else

    D = A\B;

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Grabbing the data back
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %  Extract displacements from D vector into subvectors
if Option=='B' || Option=='E'     
    TensileSlipDisp = D(1:NUM,:);         
    StrikeSlipDisp = D((NUM+1):(2*NUM),:);
    DipSlipDisp = D(((2*NUM)+1):(3*NUM),:);
elseif Option=='D'      
    StrikeSlipDisp = D(1:NUM,:);        
    DipSlipDisp = D((NUM+1):(2*NUM),:);
    TensileSlipDisp=zeros(size(StrikeSlipDisp));
elseif Option=='F'  
    TensileSlipDispE1   = D(1:NUME1,:);                          %Elastic1
    StrikeSlipDispE1    = D((NUME1+1):(2*NUME1),:);              %Elastic1
    DipSlipDispE1       = D(((2*NUME1)+1):(3*NUME1),:);          %Elastic1
    TensileSlipDispE2   = D((3*NUME1)+1:(3*NUME1)+NUME2,:);      %Elastic2 
    StrikeSlipDispE2    = D((3*NUME1)+1+NUME2:(3*NUME1)+2*NUME2,:);     %2     
    DipSlipDispE2       = D((3*NUME1)+1+(2*NUME2):end,:);        %Elastic1   
    TensileSlipDisp     =[TensileSlipDispE1;TensileSlipDispE2];
    StrikeSlipDisp      =[StrikeSlipDispE1;StrikeSlipDispE2];     
    DipSlipDisp         =[DipSlipDispE1;DipSlipDispE2];
    
elseif Option=='C' 
    %passing results into friction solver
     [TensileSlipDisp,StrikeSlipDisp,DipSlipDisp] = LinearCompFrictionSolver3d(D,A,Sf,Mu,NUM,Tnn,Tss,Tds);

end
    clear A B
                                               
if Option=='E'
    Pxx=0;
    Pyy=0;
    Pzz=0;
    Pxy=0;
    Pxz=0;
    Pyz=0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating and plotting the traction vectors.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Total traction vector
% T=sqrt((Tx.^2)+(Ty.^2)+(Tz.^2)); %Bars, see Pollard Book Eq 6.53,http://math.stackexchange.com/questions/586245/what-does-double-vertical-line-means-in-linear-algebra/586375
% %Redef normal stress
% Tn=NormalTraction;
% %Trig to find shear
% tosqrt=(T.^2)-(Tn.^2);
% TsTrig=sqrt(abs(tosqrt));
% 
% %Creating and drawing Shearing traction vector
% ShearToFill=zeros(size(FaceNormalVector));
% for i=1:size(FaceNormalVector(:,1));
% N=FaceNormalVector(i,:)';
% T=[Tx(i,:);Ty(i,:);Tz(i,:)];
% %Equations below taken in part from
% %''arrowsmith510.asu.edu/TheLectures/Lecture16/Lecture16_3Dstress.ppt''
% %Now for the shear traction; use the McKenzie construction
% B = cross(T,N); %vector normal to the plane containing T and N
% Ts = cross(N,B); %shear traction direction
% Ts_mag =  sqrt(Ts(1)^2 + Ts(2)^2 + Ts(3)^2);%TsTrig(i); 
% Ts(1) = Ts(1)./Ts_mag;
% Ts(2) = Ts(2)./Ts_mag;
% Ts(3) = Ts(3)./Ts_mag;
% ShearToFill(i,:)=Ts';
% end
% figure;quiver3(MidPoint(:,1),MidPoint(:,2),MidPoint(:,3),ShearToFill(:,1),ShearToFill(:,2),ShearToFill(:,3))
% xlabel('x'); ylabel('y'); axis('equal'); title('Shear Traction Vector');
% figure;quiver3(MidPoint(:,1),MidPoint(:,2),MidPoint(:,3),Tx(:,1),Ty(:,1),Tz(:,1))
% xlabel('x'); ylabel('y'); axis('equal'); title('Total Traction Vector');

end %end entire func

