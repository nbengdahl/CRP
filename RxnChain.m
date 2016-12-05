close all; clear all; clc;
% RxnChain.m - an example included as a supplement to "Lagrangian 
% simulation of mixing and reactions in complex geochemical systems" 
% by Engdahl, Benson, Bolster in Water Resources Research
% 
% Example run script to compute the reaction chain simulation and generate 
%  Figure 4. Please READ IT so you understand what the 
%  reaction is supposed to do.
%
% Requires:
%    1. Control m-file    (that's this file)
%       - Note: The run file "RxnChain_Auto.in" will be generated automatically 
%    2. Rxn_dBase.dat      (Database file for this example)
%  **3. PhreeqcRM library installed and added to the Matlab path
% 
%  ** You're on your own for that last one but it's not too bad.
% 
% You are free to modify and/or distribute these examples but please do
% cite the above paper. 
% 
% Direct questions to nick.engdahl@wsu.edu
% =======================================================================

% n-dimensional colocation density (needed below)
nus=@(s,dt,D1,D2,n) exp(-(abs(s).^2)/(4.*dt.*(D1 + D2)))./((4.*pi.*dt.*(D1 + D2)).^(n./2));

% >>>>>>>> DEFINE THE REACTION HERE <<<<<<<<
% Reaction: A + B => C; C + E <-> CE; E->null (kinetic, eqilibrium, decay)

kf=1.0;     % Kinetic rate constant
lam=1e-3;    % Decay rate for E

A0=1.0; % Initial concentrations
B0=1.15; 
E0=1.25; 

% ========================================================================
% This block auto-generates the required input file for PhreeqcRM 
% 
fid = fopen('RxnChain_Auto.in','w');
fprintf(fid,'%s\n','# Irreversible kinetic example');
fprintf(fid,'%s\n','SOLUTION 1');
fprintf(fid,'\t %s %f \n','[A]',A0); % Mass is only asigned to every other particle so double it
fprintf(fid,'\t %s %f \n','[E]',E0);
fprintf(fid,'%s\n','SOLUTION 2');
fprintf(fid,'\t %s %f \n','[B]',B0); 
fprintf(fid,'%s\n','KINETICS 1');
fprintf(fid,'%s\n','[AB]_rxn');
fprintf(fid,'\t %s\n','-formula [A] 1 [B] 1 [C] -1');
fprintf(fid,'\t %s %f %d %d\n','-parms',kf,1,1);
fprintf(fid,'\t %s\n','-m0 0.0');
fprintf(fid,'\t %s\n','-m 0.0');
fprintf(fid,'%s\n','KINETICS 2');
fprintf(fid,'%s\n','[AB]_rxn');
fprintf(fid,'\t %s\n','-formula [A] 1 [B] 1 [C] -1');
fprintf(fid,'\t %s %f %d %d\n','-parms',kf,1,1);
fprintf(fid,'\t %s\n','-m0 0.0');
fprintf(fid,'\t %s\n','-m 0.0');
fprintf(fid,'%s\n','[E]_decay');
fprintf(fid,'\t %s\n','-formula [E] 1.0 ');
fprintf(fid,'\t %s %f\n','-parms',lam);
fprintf(fid,'\t %s\n','-m0 0.0');
fprintf(fid,'\t %s\n','-m 0.0');
fprintf(fid,'%s\n','END');
fclose(fid);
%
% ========================================================================

if not(libisloaded('libphreeqcrm'))
    loadlibrary('libphreeqcrm','RM_interface_C.h');
    disp('Hold on, loading library...');
    pause(1);
end

% ==== Base unit conversion factor ====
c_fac=1000; % mol/L to mg/L FOR A MOLECULAR MASS OF ONE - CHANGE OTHERWISE!
% This is another place to be careful, some units are in mg/L others mol/L
% so you may need to load your masses and define species specific c_fac


% ======== DEFINE PARTICLES ========
nthreads=2;  % Number of threads to use
nprt=1500;      % Number of "cells" which are particles for me...

xL=1;
xmin=-0.5;
xp=linspace(xmin,xmin+xL,nprt);    % Make an evenly spaced line of particles
dx=xL./nprt;                        % Average distance bewteen particles

% Array for particle positions over time
Xp=zeros(nprt,2);
Xp(:,1)=xp;         % x-position
Xp(:,2)=0;          % y-position (not used yet)

% Master diffusion coefficient, here only one is used for all species but
% you could define more and call "nus" with the correct value for each
D0=1e-2;

% ======== DEFINE TIME STEPPING ========
cur_time=0.0;       % Starting time
dt=0.1;             % Time step size
t_max=1500;          % Max time

nsteps = round(t_max./dt,0);   

% ======== DEFINE PROCESSES AND BOUNDARIES ========
diffusion_on=1;     % Explicitly simulate diffusion? (1=yes)
periodic=0;         % Periodoc(1) or rigid(0) boundaries

gint=50;  % Skip interval for plotting, set to 1 to plot all times
plot_me=1;  % Enable/Disable simulation plot over time

% ======== DEFINE PHYSICAL CONDITIONS ========
% These don't make much difference for this problem
den_0=1.0;  % Density
prs_0=1.0;  % Pressure
tmp_0=25.0; % Temperature
sat_0=1.0;  % Saturation
por_0=1.0;  % Posority
vol_0=1.0;  % Volume

% Vectors for all of the above
den_v=(den_0).*ones(nprt,1);
prs_v=(prs_0).*ones(nprt,1);
tmp_v=(tmp_0).*ones(nprt,1); 
sat_v=(sat_0).*ones(nprt,1);
por_v=(por_0).*ones(nprt,1);
vol_v=(vol_0).*ones(nprt,1);

if (diffusion_on)
    D=D0./2; % If explicit diffusion, split bewteen mass transfer and Brownian motion
else
    D=D0;    % Mass trasnfer based diffusion only
end

ds=dx.*ones(nprt,1); % Vector of ds values needed below

pfid=calllib('libphreeqcrm','RM_Create',nprt,nthreads);

[LdDb]=calllib('libphreeqcrm','RM_LoadDatabase',pfid,'Rxn_dBase.dat');
if (strcmp(LdDb,'IRM_FAIL'))
    error('Database load error');
end

[~]=calllib('libphreeqcrm','RM_SetSpeciesSaveOn',pfid,1);
[~]=calllib('libphreeqcrm','RM_OpenFiles',pfid);
[~]=calllib('libphreeqcrm','RM_SetRepresentativeVolume',pfid,vol_v); % Should be 
[~]=calllib('libphreeqcrm','RM_SetSaturation',pfid,sat_v);
[~]=calllib('libphreeqcrm','RM_SetPorosity',pfid,por_v);
[~]=calllib('libphreeqcrm','RM_SetTime',pfid,cur_time);
[~]=calllib('libphreeqcrm','RM_SetTimeStep',pfid,dt);
[~]=calllib('libphreeqcrm','RM_SetDensity',pfid,den_v);
[~]=calllib('libphreeqcrm','RM_SetTemperature',pfid,tmp_v);
[~]=calllib('libphreeqcrm','RM_SetPressure',pfid,prs_v);

mask=zeros(1,nprt); mask(1)=1;
Imsk=calllib('libphreeqcrm','RM_SetPrintChemistryMask',pfid,mask);
Ipch=calllib('libphreeqcrm','RM_SetPrintChemistryOn',pfid, 0, 0, 0);

[~]=calllib('libphreeqcrm','RM_RunFile',pfid,1,1,1,'RxnChain_Auto.in');

ncomp=calllib('libphreeqcrm','RM_FindComponents',pfid);
ngrd2=calllib('libphreeqcrm','RM_GetGridCellCount',pfid);
if (ngrd2~=nprt)
    error('Cell/particle count changed');
end
nspec=calllib('libphreeqcrm','RM_GetSpeciesCount',pfid);
nchem=calllib('libphreeqcrm','RM_GetChemistryCellCount',pfid);


% =======================================================================
%                 INITIAL CONDITION
% =======================================================================
% === % (1) SOLUTIONS, (2) EQUILIBRIUM_PHASES, (3) EXCHANGE, 
% === % (4) SURFACE, (5) GAS_PHASE, (6) SOLID_SOLUTIONS, and (7) KINETICS
 
iniC1=(-1).*ones(nprt,7); 
iniC1(:,7)=1; 

rxn_frac=0.2;   % How abundant is the E decay reaction (this sets 20%)
n_rxn=round(nprt.*rxn_frac,0);
if (n_rxn<=0)
    error('No sites for reaction, increase fraction');
end

rng(317821);
rn=randperm(nprt);
iniC1(rn(1:n_rxn),7)=2;

% Split (mixing limited) initial condition
iniC1(:,1)=1.*(xp<0) + 2.*(xp>=0);  % Half get Sln 1 others Sln 2

% % or use the well-mixed initial comdition
% iniC1(:,1)=1;
% iniC1(2:2:nprt,1)=2;
%
%            >>>>> End of initial condition <<<<<
% =======================================================================

% Send initial conditions/definitions to the solver
[IOin]=calllib('libphreeqcrm','RM_InitialPhreeqc2Module',pfid,iniC1,[],[]);

if (1==1) % Warning: This CAN crash Matlab if there is a file error...
% Make the lists of the components/species 

species_list=cell(nspec,1);
s_name=char('string');
for n=0:nspec-1 % Index starts at zero
[stat2,s_name]=calllib('libphreeqcrm','RM_GetSpeciesName',pfid,n,s_name,10);
species_list{n+1}.name=s_name;
end
compnts_list=cell(ncomp,1);
for n=0:ncomp-1 % Index starts at zero
[stat2,s_name]=calllib('libphreeqcrm','RM_GetComponent',pfid,n,s_name,10);
compnts_list{n+1}.name=s_name;
end
compnts_list_s={' '};
species_list_s={' '};
compnts_list_s=repmat(compnts_list_s,ncomp,1);
species_list_s=repmat(species_list_s,nspec,1);
for m=1:ncomp
compnts_list_s{m}=char(compnts_list{m}.name);
end
for m=1:nspec
species_list_s{m}=char(species_list{m}.name);
end

end
% return        % Break point for checking which species have been loaded 

c=zeros(nprt,ncomp);
species_c=zeros(nprt,nspec);
[~,c]=calllib('libphreeqcrm','RM_GetConcentrations',pfid,c);
[~,species_c]=calllib('libphreeqcrm','RM_GetSpeciesConcentrations',pfid,species_c);

% Store initial concentrations
species_c0=species_c;
c0=c;

% ========================================================================
%               Take the first mass transfer step now
%
%  -- NOTE: For large particle counts it is more efficient to use a search
%  tree and only evaluate the particles that are close enough to colocate;
%  however, the method used below doesn't require any external functions.
%
cur_time = cur_time + dt;
[~]=calllib('libphreeqcrm','RM_SetTime',pfid, cur_time);  %   Update current time

X=repmat(Xp(:,1)',nprt,1); % Matrix of all x positions
Y=repmat(Xp(:,2)',nprt,1); % Matrix of all y positions
S=sqrt(((X-transpose(X)).^2) + ((Y-transpose(Y)).^2));  % Separation dist.
if (periodic)
S=S.*(S<=(xL./2)) + (S-xL).*(S>(xL./2));
end
Nu=nus(S,dt,D,D,1);      % Colocation density (n-d)
for m=1:nspec
inC=species_c(:,m);
M=repmat(species_c(:,m)',nprt,1); % Matrix of the masses
Md=M-transpose(M);       % Difference in mass
dM=((0.5).*(Md.*Nu)); 
ddM=dM*ds;
species_c(:,m)=species_c(:,m)+ddM;
end
[~]=calllib('libphreeqcrm','RM_SpeciesConcentrations2Module',pfid, species_c);  % Send concentrations 

% Execute the reaction step
[IOt]=calllib('libphreeqcrm','RM_RunCells',pfid);

[~,c]=calllib('libphreeqcrm','RM_GetConcentrations',pfid,c);
[~,species_c]=calllib('libphreeqcrm','RM_GetSpeciesConcentrations',pfid,species_c);
%
%                ----- END of initial step -----
% ========================================================================

% Place holders for output concentrations
c2=c.*0;    
species_c2=species_c.*0; 

s_lst=[6 7 8 9 10];  % List to define which species to save, make 1:nspec for all

if (plot_me)
figure;
title('Species concentrations');
for mm=s_lst
plot(Xp(:,1),species_c(:,mm).*c_fac,'.'); hold on; grid on;
end
legend(species_list_s(s_lst)); xlabel('x-position');ylabel('Concentration');
drawnow; pause(0.25);
end

% Output storage and initial concentrations
t=dt.*(1:nsteps);
C_t=zeros(nsteps,length(s_lst));    % Species over time
C_t0=zeros(1,length(s_lst));        % Initial species concentrations

% Log the first reaction step
nn=1;
for jj=s_lst
    % --> c_fac NEEDS TO BE MODIFIED FOR DIFFERENT REACTIONS/DATABASES <--
    C_t(1,nn)=sum(species_c(:,jj).*c_fac.*ds); % Convert to mg/L (assumes all species have unit mass)
    nn=nn+1;
end

% Loop over the remaining steps
for isteps=2:nsteps
    
% =======================================================================    
% Take a diffusive step and enforce periodic boundaries
if (diffusion_on)
rd=randn(nprt,1);
Xp(:,1)=Xp(:,1) + sqrt(2.*(D).*dt).*rd;

if (periodic)
    % Logical implementation of the periodic boundary
    Xp(:,1)=Xp(:,1).*((Xp(:,1)<=max(xp))&(Xp(:,1)>=min(xp))) + ...
    (max(xp)+(min(xp)-Xp(:,1))).*(Xp(:,1)<=min(xp)) + ...  
    (min(xp)+(Xp(:,1)-max(xp))).*(Xp(:,1)>=max(xp));
else
    % Logical implementation of the bounce back boundary
    Xp(:,1)=Xp(:,1).*((Xp(:,1)<=max(xp))&(Xp(:,1)>=min(xp))) + ...
    (min(xp)+(min(xp)-Xp(:,1))).*(Xp(:,1)<=min(xp)) + ...  
    (max(xp)+(max(xp)-Xp(:,1))).*(Xp(:,1)>=max(xp));
end
end

% ======== Mass exchange step =========
X=repmat(Xp(:,1)',nprt,1); % Matrix of X positions
Y=repmat(Xp(:,2)',nprt,1); % Matrix of Y positions
S=sqrt(((X-transpose(X)).^2) + ((Y-transpose(Y)).^2));  % Separation dist.
if (periodic)
S=S.*(S<=(xL./2)) + (S-xL).*(S>(xL./2));
end
Nu=nus(S,dt,D,D,1);      % Colocation density (n-d)

for m=1:nspec
inC=species_c(:,m);
M=repmat(species_c(:,m)',nprt,1); % Matrix of the masses
Md=M-transpose(M);       % Difference in mass
dM=((0.5).*(Md.*Nu));
ddM=dM*ds;
species_c(:,m)=species_c(:,m)+ddM;
end

% =======================================================================
% *** REACTION STEP
% 
[~]=calllib('libphreeqcrm','RM_SpeciesConcentrations2Module',pfid, species_c);      %  Concentrations
cur_time = cur_time + dt;
[~]=calllib('libphreeqcrm','RM_SetTime',pfid, cur_time);     %  Current time

% Run cells with the now transported concentrations
[IOr]=calllib('libphreeqcrm','RM_RunCells',pfid);  

% Transfer data back from PhreeqcRM for transport
[~,c2]=calllib('libphreeqcrm','RM_GetConcentrations',pfid, c2);
[~,species_c2]=calllib('libphreeqcrm','RM_GetSpeciesConcentrations',pfid, species_c2); 
[~,den_v]=calllib('libphreeqcrm','RM_GetDensity',pfid, den_v);
[~,vol_v]=calllib('libphreeqcrm','RM_GetSolutionVolume',pfid, vol_v); 	
c=c2;
species_c=species_c2;
% 
% *** END OF REACTION STEP
% =======================================================================

if (plot_me)&& (ceil(isteps/gint)==(isteps/gint))
    cla;
    title(['Species concentrations at t=',num2str(cur_time)]);
    for mm=s_lst
    plot(Xp(:,1),species_c(:,mm).*c_fac,'.'); hold on; grid on;
    end
    plot([xmin xmin+xL],(A0).*[1 1],'--m'); ylim([0 A0]);
    legend(species_list_s(s_lst)); xlabel('x-position');ylabel('Concentration');
    drawnow; pause(0.25);
end

% Log the species concentrations for this time
nn=1;
for jj=s_lst
    C_t(isteps,nn)=sum(species_c(:,jj).*c_fac.*ds);
    nn=nn+1;
end

end

%%
% close all;

figure;  
t2=logspace(-2,2,100);
for jj=1:length(s_lst) % Just plot A
plot(t,C_t(:,jj),'o','MarkerSize',3); hold on; grid on;
end

set(gca,'xscale','log','yscale','log'); xlim([1e-2 1e3]); ylim([1e-3 1e0]);
ylim([A0./1000 A0]);
xlabel('Time');ylabel('Concentration');
legend(species_list_s(s_lst),'location','NorthWest');

% xlim([1e-1 1e3]); ylim([1e-3 1e0]);

%%
% Unload the library
[IDEL]=calllib('libphreeqcrm','RM_Destroy',pfid);

unloadlibrary('libphreeqcrm');