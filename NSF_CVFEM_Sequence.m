%% NSF_CVFEM_Sequence.m
%Modified to work with Sequence Outputs

%% Front Matter
%Simulated advancing simplified clinoform-Field Scale

tic
%This code was developed as part of a US National Science Foundation Project  
%Collaborative Research: Exploring the linkages between Sea-Level Change, Sediment Transport and %Geomorphology on Coastal Freshwater Water Sequestration
%Grant Number NSF 1925506

%The principal coder was Vaughan R. Voller, University of Minnesota, volle001@umn.edu
%With contributions from Mark Person, Nafis Sazeed, and Loc Luong, New Mexico Tech. 
%The particular code has been modified by Nafis Sazeed to work with Sequence outputs with significant contributions from Vaughan R. Voller, Mark Person and Loc Luong.

%This version has been prepared as supplementary material for the paper:-  

%The effects of late pleistocene sea-level fluctuations and sediment transport processes on the sequestration of fresh and brackish water in continental shelf environments

% Nafis Sazeed, Mark A. Person, Vaughan R. Voller, Michael S. Steckler, Loc Luong, Eric Hutton, Kerry Key, Huy Le, Celine Jo Grall

%Submitted to Water Resources Research, 2025

%% Code starts

clear all
close all

%% Importing Sequence Output File
ncfile='case54.nc';
%Run the following code to import the sequence output file NS
sequenceImport

%% Basic Settings
%Basic Geometric data--set for base-case in wrr paper June 2022
Breath=seq_nodes*(x_of_cell(1,2)-x_of_cell(1,1));  %Breath (m) of domain (Length of long profile)
cols=seq_nodes-2; % Number of grid columns--equally spaced
Delx=Breath/(cols-1); %fixed width between columns
xcols=linspace(0,Breath,cols); %xlocations of columns

% delt=100; % Time step days
delt=1000*365; %1000 Years

%Physical properties
consea=1;  %Salt concentration in sea
rhorel=1.025;%Relative density of saturated saline rhosat/rhowater

eps=0.35; %Aquifer porosity
stov=0.001; %0.001; %Aquifer storage

anis=100; %kx/kz

%Hydraulic conductivity m/day, othotropic
kxval=0.0086/1000;
kyval=kxval/anis;

ksf=8.6;%8.6*0.2;
kg=8.6*0.8;%ksf*0.8;
kse=kxval*20;

%diffusion and dispersion coefficients to be used in anisotropic dispersion
Dmol= 0.00003;%  molecular diffusion  m^2/day
aL=50; % Longitudinal (m)
aT=aL/10; %Transverse

%% SEQUENCE IMPORT
totstep=tot_tstep;
sealevelin=sea_level(1,1);%initial sealevel from sequence

%% INTERPOLATING DIMENSION by NS
%We will interpolate the top surface, bedrock and sediment surface
%to a fixed number of columns, so that we can use the same grid for all cases
new_columns=200;
xq=linspace(0,Breath,new_columns);

topsurface(1,:)=topsurface(2,:);

for i=1:totstep
    vq1=interp1(xcols,topsurface(i,:),xq);
    vq2=interp1(xcols,bedrock_position(i,:),xq);
    vq3=interp1(xcols,sed_surf(i,:),xq);
    eta(i,:)=vq1;
    etabot(i,:)=vq2;
    etabed(i,:)=vq3;
end
cols=max(size(xq));
Delx=Breath/(cols-1);
xcols=xq; %locations of columns 
Dely=2; %vertical distance between nodes in each column 

%% Grid Initial domain

ytop=etabed(1,:);  %initial elevation of top domain nodes
ybot=etabed(1,:)-4*Dely; %initial elevation of bottom domain nodes

ncol=round((ytop-ybot)/Dely)+1; %number of nodes in each column
N=sum(round((ytop-ybot)/Dely)+1); %number of nodes in domian

%storage and initial values for head values
phi=max(ytop)*ones(N,1);
phi_insert=zeros(1,cols);

%storage and initial values for concentration
con_insert=zeros(1,cols);
age_insert=zeros(1,cols);

%set salt intitial cons of points vertically below sea level to value 1.
con=[];
age_=[];
for jj=1:cols
    if ytop(jj)>sealevelin
        con=[con;zeros(ncol(jj),1)];
        age_=[age_;zeros(ncol(jj),1)];
    else
        con=[con;ones(ncol(jj),1)];
        age_=[age_;zeros(ncol(jj),1)];
    end
end

%% Time Stepping

tolh=1e-5; %convergence tolerance for head
tols=1e-7; %for concentration 
tola=1e-4; %for age

printval=1; %50;   %paramters to control plotting through time
printspace=1; %50;

for tstep=1:500 %Main time loop
    
    %% Adjustment of strat and set sea-level
    ytopold=ytop;
    ybotold=ybot;
    if tstep<totstep+1
        ytop=eta(tstep,:);
        ybot=etabed(tstep,:)-4*Dely;
    end
    
    %Sometimes Sequence output has drastic changes in topography at the first timestep. The following line gets rid of this issue through hardwiring
    if (tstep==1); 
        ytop=etabed(tstep,:);
    end

    %We use the code plotStratigraphy to plot the stratigraphy and import the shapes as x and y values with the suffix se(shelfedge), sf(shoreface), and g(fluvial)
    pl_tstep=tstep+1;
    plotStratigraphy
    
    ncolpre=ncol; %store previous value of nodes on each column 
    ncol=round((ytop-ybot)/Dely)+1; %calculate new value. Will only change if
    %top nop node is at a distance > 1.5 Dely
    % above its neigbour below
    Btop=[]; %list of nodes on top surface
    Btopsea=[];  %list of top domain nodes under sea
    x=[]; %list of x locations of domain nodes
    y=[]; %list of y locations of domain nodes
    for jj=1:cols
        x=[x;Delx*(jj-1)*ones(ncol(jj),1)];
        y=[y;linspace(ybot(jj),ybot(jj)+(ncol(jj)-2)*Dely,ncol(jj)-1)';ytop(jj)];
        %Find Btop, Btopsea Nodes on top surface, top surface nodes under sea
        Btop=[Btop,size(y,1)];
        if ytop(jj)<sea_level(tstep,1) %XXXX
            Btopsea=[Btopsea,size(y,1)];
        end
    end
    Bbot=[1,Btop(1:cols-1)+1]; %List of bottom nodes
    
    N=size(x,1); %updated number of node points
    
    % Putting in the d_eta term by NS
    del_eta=(ytop-ybot)-(ytopold-ybotold); %Calculate changes in elevation for all columns at each timestep
    d_eta=[];
    for i=1:cols
        for j=1:ncol(i)
            d_eta=[d_eta;del_eta(i)];
        end
    end
    
    %store values from previous grid  XXXXXX
    phipre=phi;
    conpre=con;
    agepre=age_;
    
    
    %resize phi and con as grid expands
    con=[];
    phi=[];
    age_=[];
    last=0;
    for jj=1:cols
        first=last+1;
        last=first+ncolpre(jj)-1;
        
        
        %No change iin number of nodes in columne
        if ncolpre(jj)==ncol(jj)
            phi=[phi;phipre(first:last)];
            con=[con;conpre(first:last)];
            age_=[age_;agepre(first:last)];
        end
        
        %Remove node in columne due to errosion
        if ncolpre(jj)>ncol(jj)
            ndiff=ncolpre(jj)-ncol(jj);
            % Using ndiff instead of 1 node to erode multiple nodes
            
            phi=[phi;phipre(first:last-ndiff)];
            con=[con;conpre(first:last-ndiff)];
            age_=[age_;agepre(first:last-ndiff)];
        end
        
        if ncolpre(jj)<ncol(jj)  %This condition indicates that
            % node points have been added
            
            %CODE TO ADD MULTIPLE NODES
            ndiff=ncol(jj)-ncolpre(jj);
            
            %       ADDED The n_diff code
            %assign current values to nodes up to lowest inserted node
            phi=[phi;phipre(first:last-1)];
            con=[con;conpre(first:last-1)];
            age_=[age_;agepre(first:last-1)];
            
            %We then overwrite the added node value---
            %using the stored interpolatd values, calculated immediately after
            %the last time step solution
            mnode=size(con,1); % gives top node position in current column
            dy=y(mnode+ndiff+1)-y(mnode);
            
            %interpolation ratio
            phirat=(phipre(last)-phipre(last-1))/dy;
            conrat=(conpre(last)-conpre(last-1))/dy;
            agerat=(agepre(last)-agepre(last-1))/dy;
            
            for ik=1:ndiff
                phi=[phi;phipre(last-1)+phirat*Dely*ik];
                con=[con;conpre(last-1)+conrat*Dely*ik];
                age_=[age_;agepre(last-1)+agerat*Dely*ik];
            end
            
            phi=[phi;phinew(last)];
            con=[con;connew(last)];
            age_=[age_;agenew(last)];
            
        end
        
    end
   
    %% Make Mesh
    
    %We bring in two adjacent columns at a time and then use the
    %inbuilt MATLAB routine delaunay to make the elemnts in the mesh
    t=[];
    first=1;
    for jj=1:cols-1
        last=first+ncol(jj)+ncol(jj+1)-1;
        xcc=x(first:last);
        ycc=y(first:last);
        tcc=delaunay(xcc,ycc)+first-1;
        t=[t;tcc];
        first=first+ncol(jj);
    end
    %The resulting array 't' has one line for each element in mesh.
    %this line contains a list of the three node inecies
    %(in couter clockwise order) that are at the vertices on the element
    
    %NOTE the MATALB command
    %triplot(t,x,y, '-b');
    %will print out the grid at any time asked for
    
    Ntri=size(t,1);  %the size of t is the number of triangle elements
    
    %% Geometric properties of mesh
    %Here we calculate some geometric properties of the mesh
    %And asseemble the support array sup
    
    Volp=zeros(N,1);    %CV volume
    xmid=zeros(Ntri,1); %elemnt mid point
    ymid=zeros(Ntri,1);
    Nx=zeros(Ntri,3);   %Derivatives of shape functions
    Ny=zeros(Ntri,3);
    
    for itri=1:Ntri
        k1=t(itri,1); %global number of 1st node in trinagle itri
        k2=t(itri,2); %2nd node
        k3=t(itri,3); %3rd node
        %element volume
        v=(x(k2)*y(k3)-x(k3)*y(k2)-x(k1)*y(k3)+x(k1)*y(k2) ...
            +y(k1)*x(k3)-y(k1)*x(k2))/2;
        
        
        Volp(k1)=Volp(k1)+v/3; %contribution to control volume
        Volp(k2)=Volp(k2)+v/3;
        Volp(k3)=Volp(k3)+v/3;
        
        xmid(itri)=(x(k1)+x(k2)+x(k3))/3; %mid point of elemnt
        ymid(itri)=(y(k1)+y(k2)+y(k3))/3;
        
        %derivatives of shape functions
        Nx(itri,1)= (y(k2)-y(k3))/(2*v);   %Nx=shape function derivative
        Nx(itri,2)= (y(k3)-y(k1))/(2*v);   %the index 1, 2 or 3
        Nx(itri,3)= (y(k1)-y(k2))/(2*v);   %refers to local tri elemnt node
        Ny(itri,1)=-(x(k2)-x(k3))/(2*v);   %Ny=shape function derivative
        Ny(itri,2)=-(x(k3)-x(k1))/(2*v);
        Ny(itri,3)=-(x(k1)-x(k2))/(2*v);
        
    end
    
    if min(xg)>0
        xg(xg==min(xg))=0;
    end
    
    if max(yg)<max(y)
        yg(yg==max(yg))=max(y);
    end
    
    %We will now create polygonal boundaries for the different facies from Sequence- by NS
    yqq=[ybot+4*Dely,ybot(cols:-1:1)];
    xqq=[xq,xq(cols:-1:1)];
    
    in_mat_se=inpolygon(xmid,ymid,xse*1000,yse);
    in_mat_sf=inpolygon(xmid,ymid,xsf*1000,ysf);
    in_mat_g=inpolygon(xmid,ymid,xg*1000,yg);
    % in_mat_bed=inpolygon(x,y,xqq,yqq);
    
    %% SET BOUNDARY CONDITIONS
    
    %BOUNDARY NODE POINTS
    
    %We are assumming that only no flow or prescibed bounday conditions
    %are applied.
    
    %By the manner in which our CVFEM discrete systesm is constructed no-flow
    %boundaries are the defualy setting.
    
    % We set precribed boundary conditiosn as follows
    % The form of our discrte equation for a given node is (eq 12)
    
    % (mu*V + delt*a_i)*u_i =mu*V a^old_i + delt*sum(a_j*u_j)
    
    % To account for boundaries we add additionl coeeficents and terms
    
    % (mu*V + delt*a_i +BCu_i)*u_i =mu*V a^old_i + delt*sum(a_j*u_j) +BBu_i
    
    %The one-d arrays BCu and BBu are initialized with zeros
    %at nodes where a fixed value 'val' needs to be applied we set
    %BCu_i=1e18
    %BBu_i=1e18*val
    % ensuring, on solution, of the discrete equation that the correct value
    % is imposed at node i
    
    
    %Btop--stored node numbers on top boundary
    %Btopsea--stored top nodes at or below sea level.
    
    %Arrays for Boundary conditions
    %Recall
    %Btop--stored node numbers on top boundary
    %Btopsea--stored top nodes at or below sea level.
    
    Big=1e18;
    BCh=zeros(N,1); %boundary coeff values for head and solute
    BCs=zeros(N,1);
    BBh=zeros(N,1); %fixed head values
    BBs=zeros(N,1); %fixed solute consentration values
    
    %Fixed head value on top boundary
    BCh(Btop,1)=Big;
    BBh(Btop,1)=Big*(y(Btop));
    %correctiion for sea water
    BBh(Btopsea,1)=Big*(y(Btopsea)-(y(Btopsea)-sea_level(tstep,1))*(rhorel));
    
    %concentration Fixed values 0 above sealevel 1 below
    BCs(Btop,1)=Big;
    BBs(Btopsea)=consea*Big;
    %Note for the submerged nodes stored in Btopsea we will
    %overwrite this condition when we detect outflow--see below
    
    %age boundary conditions
    BBa=zeros(N,1);
    BCa=zeros(N,1);
    BBa(Btopsea,1)=Big;
    BCa(Btop,1)=Big;
    BBa(Btopsea)=Big;
    
    %% Head Solution
    
    %coefficents
    BBvar=zeros(N,1); %Variable density source
    
    kx=ones(Ntri,1)*kxval; %x direction conductivity value
    ky=ones(Ntri,1)*kyval; %y direction conductivity value
    
    AA=sparse(N,N);
    BB=zeros(N,1);
    %VV
    
    
    for itri=1:Ntri
        % Conductivity values for each element according to their facies type
        if(in_mat_se(itri)==1)
            kx(itri)=kse;
            ky(itri)=kse/anis;
        elseif(in_mat_sf(itri)==1)
            kx(itri)=ksf;
            ky(itri)=ksf/anis;
        elseif(in_mat_g(itri)==1)
            kx(itri)=kg;
            ky(itri)=kg/anis;
        end
        
        
        %The element has three vertices nodes
        %Given local indices 1, 2 and 3
        %We will call at each node in the element in turn
        %and build the coefficents assocaited with that node
        %when we work with node 1 the neigbours are 2 and 3 (counter-clock)
        %when we work with node 2 the neigbours are 3 and 1 (counter-clock)
        %when we work with node 3 the neigbours are 1 and 2 (counter-clock)
        %We store this ordering in
        cyc=[1,2,3;2,3,1;3,1,2];
        
        for node=1:3
            ii=cyc(node,1);  %local node
            jj=cyc(node,2);
            kk=cyc(node,3);
            
            k1=t(itri,ii); %global node number
            k2=t(itri,jj);
            k3=t(itri,kk);
            
            Nx1=Nx(itri,ii);   %Nx=shape function derivative for elemnt
            Nx2=Nx(itri,jj);
            Nx3=Nx(itri,kk);
            Ny1=Ny(itri,ii);   %Ny=shape function derivative
            Ny2=Ny(itri,jj);
            Ny3=Ny(itri,kk);
            
            %Face1
            delx= (x(k1)+x(k2)+x(k3))/3-(x(k1)+x(k2))/2;
            dely= (y(k1)+y(k2)+y(k3))/3-(y(k1)+y(k2))/2;
            
            face1_k1=kx(itri)*Nx1*dely-ky(itri)*Ny1*delx;
            face1_k2=kx(itri)*Nx2*dely-ky(itri)*Ny2*delx;
            face1_k3=kx(itri)*Nx3*dely-ky(itri)*Ny3*delx;
            
            %varibale density source(face value)
            BBvar(k1)=BBvar(k1)-ky(itri)*(((rhorel-1))/12)*(5*con(k1)+5*con(k2)+2*con(k3))*delx;
            
            
            %Face2
            delx= -(x(k1)+x(k2)+x(k3))/3+(x(k1)+x(k3))/2;
            dely= -(y(k1)+y(k2)+y(k3))/3+(y(k1)+y(k3))/2;
            
            
            face2_k1=kx(itri)*Nx1*dely-ky(itri)*Ny1*delx;
            face2_k2=kx(itri)*Nx2*dely-ky(itri)*Ny2*delx;
            face2_k3=kx(itri)*Nx3*dely-ky(itri)*Ny3*delx;
            
            %Ax=b Solver
            AA(k1,k1)=AA(k1,k1)+delt*(face1_k1+face2_k1);
            AA(k1,k2)=AA(k1,k2)+delt*(face1_k2+face2_k2);
            AA(k1,k3)=AA(k1,k3)+delt*(face1_k3+face2_k3);
            %XXXXX
            
            %varibale density source(face value)
            BBvar(k1)=BBvar(k1)-ky(itri)*((rhorel-1)/12)*(5*con(k1)+2*con(k2)+5*con(k3))*delx;
        end
    end
    
    
    
    %Update Head
    
    %Uses Jacobi (works well with vectorization)
    %provides solution of equation 12 in text--written as
    %(mu*V + delt*a_i +BCu_i)*u_i =mu*V a^old_i + delt*sum(a_j*u_j) +BBu_i
    %NOTE RHS=sum(asup.*phinew(sup),2)=sum(a_j*u_j)
    %and addition of varibale density souce
    %(solution is for freshwater head)
    
    % Storage set to 0 at first time step
    sto=stov;
    if tstep==1
        sto=0;
    end
    
    %VV
    Badd= -BCh-sto*Volp;
    AA = spdiags(spdiags(AA,0)+Badd,0,AA);
    BB=-sto*Volp.*phi-delt*BBvar-BBh-d_eta.*Volp*1.3*sto;
    
    phinew= AA\BB; %sparse solver
    
    %calculate volume of flow stored per time step
    %need this term to account for storage in solute transport
    Vstore=sto*(Volp.*(phinew-phi));
    
    phi=phinew; %old for new
    
    
    %%   Solute solution
    
    Netflux=zeros(N,1);
    
    AA=sparse(N,N);
    BB=zeros(N,1);
    %VV
    
    for itri=1:Ntri
        
        %Diffusion
        cyc=[1,2,3;2,3,1;3,1,2];
        for node=1:3
            ii=cyc(node,1);
            jj=cyc(node,2);
            kk=cyc(node,3);
            
            k1=t(itri,ii); %global node number of element vertices
            k2=t(itri,jj);
            k3=t(itri,kk);
            
            Nx1=Nx(itri,ii);   %Nx=shape function derivative
            Nx2=Nx(itri,jj);
            Nx3=Nx(itri,kk);
            Ny1=Ny(itri,ii);   %Ny=shape function derivative
            Ny2=Ny(itri,jj);
            Ny3=Ny(itri,kk);
            
            
            if node==1 %Dispersion for element
                
                % --contribution to discharge from fresh-water head
                qxval=-kx(itri)*(Nx1*phi(k1)+Nx2*phi(k2)+Nx3*phi(k3));
                qyval=-ky(itri)*(Ny1*phi(k1)+Ny2*phi(k2)+Ny3*phi(k3));
                
                %actual discharges at element midpoints
                %used to calculate dispersion tensor
                %also asummes rel density contribution is constant in element
                
                qxmid=qxval;
                qymid=qyval-ky(itri)*((rhorel-1)/3)*(con(k1)+con(k2)+con(k3));
                
                qx(itri)=qxmid;
                qy(itri)=qymid;
                
                qx2=qxmid^2;
                qy2=qymid^2;
                qabs=sqrt(qx2+qy2);
                
                %Dispersion from Bear 1972 SAME AS above
                Dxx=aT*qabs+(aL-aT)*qx2/qabs+Dmol*eps;
                Dyy=aT*qabs+(aL-aT)*qy2/qabs+Dmol*eps;
                Dxy=(aL-aT)*(qxmid)*(qymid)/qabs;
            end
            
            %Face1
            %dischage across face 1--Note qyface included rel. density
            qxface=qxval;
            qyface=qyval-ky(itri)*((rhorel-1)/12)*(5*con(k1)+5*con(k2)+2*con(k3));
            delx= (x(k1)+x(k2)+x(k3))/3-(x(k1)+x(k2))/2;
            dely= (y(k1)+y(k2)+y(k3))/3-(y(k1)+y(k2))/2;
            
            %contribution to flux across face due to dispersion
            face1_k1=(Dxx*Nx1+Dxy*Ny1)*dely-(Dyy*Ny1+Dxy*Nx1)*delx;
            face1_k2=(Dxx*Nx2+Dxy*Ny2)*dely-(Dyy*Ny2+Dxy*Nx2)*delx;
            face1_k3=(Dxx*Nx3+Dxy*Ny3)*dely-(Dyy*Ny3+Dxy*Nx3)*delx;
            
            %contribution to flux across fase due to flow
            %uses upwind, flow caries vatible in upstram direction
            qout=qxface*dely-qyface*delx; %flow out of vol k
            if qout>=0
                face1_k1=face1_k1-qout;
            else
                face1_k2=face1_k2-qout;
            end
            
            Netflux(k1) = Netflux(k1) + qout;
            
            %Face2
            qxface=qxval;
            qyface=qyval-ky(itri)*((rhorel-1)/12)*(5*con(k1)+2*con(k2)+5*con(k3));
            
            delx= -(x(k1)+x(k2)+x(k3))/3+(x(k1)+x(k3))/2;
            dely= -(y(k1)+y(k2)+y(k3))/3+(y(k1)+y(k3))/2;
            
            face2_k1=(Dxx*Nx1+Dxy*Ny1)*dely-(Dyy*Ny1+Dxy*Nx1)*delx;
            face2_k2=(Dxx*Nx2+Dxy*Ny2)*dely-(Dyy*Ny2+Dxy*Nx2)*delx;
            face2_k3=(Dxx*Nx3+Dxy*Ny3)*dely-(Dyy*Ny3+Dxy*Nx3)*delx;
            
            %upwind
            qout=qxface*dely-qyface*delx; %flow out of vol k
            if qout>=0
                face2_k1=face2_k1-qout;
            else
                face2_k3=face2_k3-qout;
            end
            
            Netflux(k1) = Netflux(k1) + qout;
            
            %Ax=b Solver
            AA(k1,k1)=AA(k1,k1)+delt*(face1_k1+face2_k1);
            AA(k1,k2)=AA(k1,k2)+delt*(face1_k2+face2_k2);
            AA(k1,k3)=AA(k1,k3)+delt*(face1_k3+face2_k3);
            
        end
        
    end
    
    %XXXXXXXXXXXXXX
    %Alternate Boundary Condition for out flows
    %As a default we set the concetration at nodes on the
    %seafloor and forset to  a prescribed consatnt value consea
    %this is achived with the setting
    BBs(Btopsea)=Big*consea;
    
    %But when there is outflow at a node on the seafloor,
    %we need to reset the condition as follows:
    %The nodal index of seafloor nodes are stored in
    % Btopsea
    %Thus the total number of seafloor nodes is given by
    %     nsub=size(Btopsea,2);
    %    %We can sweep over these nodes  with the loop
    %     for ii=1:nsub
    %         %Obtaining for each ii the node index of the seafloor node
    %         bnode=Btopsea(ii);
    %         %Rogniziing the column storage nature of our grid
    %         %we known that the index of the node immediatley below is bnode-1
    %         %Thus we can regonize outflow by comapring the head values
    %         % at bnode and bnode -1, i.e.,
    % %         if phi(bnode-1,1)-phi(bnode,1)>0 % we can expect outflow
    %             % and reset our boundary condition as
    % %             BBs(bnode)=Big*con(bnode-1);
    %             %i.e., with at outflow node we set the bondary concentration
    %             %equal to the value of the node immediatley below
    %             %Note this setting is expclit, lagging behind one time step
    %
    %             qvertsig=-(phi(bnode,1)-phi(bnode-1,1))/(y(bnode,1)-y(bnode-1,1))...
    %            -con(bnode-1)*(0.025);  %Full flow out check
    % %         qvertsig=-(phi(bnode,1)-phi(bnode-1,1)); %alt flow out check
    %         if qvertsig>0 % we will have outflow and thus we set
    %             BBs(bnode)=Big*con(bnode-1);
    %         end
    %     end
    
    %VV
    Badd= -BCs-eps*Volp-Vstore;
    
    AAcon = spdiags(spdiags(AA,0)+Badd,0,AA);
    BBcon=-eps*Volp.*con-BBs-d_eta.*Volp*1.3*sto.*con;
    connew= AAcon\BBcon; %sparse solver
    %connew=min(connew,1); %XXXX mass Limiter If needed
    con=connew;
    
    BCa(Btop,1)=Big;
    BBa(Btopsea,1)=0;
    
    Badd= -BCa-eps*Volp-Vstore;
    AAage= spdiags(spdiags(AA,0)+Badd,0,AA);
    BBage=-eps*Volp.*age_-BBa-delt.*Volp*eps-d_eta.*Volp*1.3*sto.*age_;
    
    agenew=AAage\BBage;
    age_ = agenew;
    
   tstep %Printing time step number
    
    % Saving output every 25th time step
    if mod(tstep,25)==0
        filename_x=sprintf("Case54hc_%d",tstep);
        save(filename_x);
    end
    

end
toc
