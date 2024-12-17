function [tsteps,nsol] = TRANSPORT2d (dt0,Tend,ifshow)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DATA
%% set-up or load load fluxes/velocities
%vx = load('vx.dat');[dum,ny] = size(vx);
%GRID hard-coded
nx = 100; 
ny = 100;
%
x0 = -1; xend = 1; 
y0 = -1; yend = 1;
z0 = 0; zend = 1;
%%
[dx,x,xc,nx,hx] = setup_grid1d (nx,x0,xend);
nz = 1;
[dy,y,yc,ny,hy] = setup_grid1d (ny,y0,yend);
[dz,z,zc,nz,hz] = setup_grid1d (nz,z0,zend);

%DATA set up velocities for testing
vx = zeros(nx+1,ny,nz)+1;
for j=1:nx+1
    for k=1:ny
        vx(j,k,1) = -1;
        %vx(j,k,1) = y(k);  
    end
end

%%make it interesting
hn1 = floor(nx/3*2); 
hn2 = floor(nx/5*4);
hny = floor(ny/2);
vx(hn1:hn2,hny,1) = 1.05;
%end interesting

vy = zeros(nx,ny+1,nz)+1; 
for j=1:nx
    for k=1:ny+1
        vy(j,k,1)= -1; 
        %vy(j,k,1)= -x(j); 
    end
end

vz = zeros(nx,ny,nz+1);

[xxc,yyc] = meshgrid(xc,yc); %% cell centers for plotting in 2d
xxc = xc; yyc = yc;

%VARIABLES: change velocity to flux

for j=1:ny
    for k=1:nz
        vx(:,j,k) = vx(:,j,k)*dy(j)*dz(k);
    end
end

for i=1:nx
    for k=1:nz
        vy(i,:,k) = vy(i,:,k)*dx(i)*dz(k);
    end
end

%VARIABLES
%% calculate (uvx,uvy) for CFL check
uvx = zeros(nx,ny); 
uvy = zeros(nx,ny);

for l=1:nx
    uvx(l,:) = 1/2*(vx(l,:,:)+vx(l+1,:,:));
end

for m=1:ny
    uvy(:,m) = 1/2*(vy(:,m,:)+vy(:,m+1,:));
end 

%umag = abs(uvx);
umag = sqrt(uvx.^2+uvy.^2);
umax = max(umag(1:end));

%VISUALIZE flow
if ifshow
%plot(dx,umag);
contourf(xxc,yyc,umag');
end

%% load ccof; should be positive
%DATA
void = ones(nx,ny,nz);
ccof = ones(nx,ny,nz); % load('ccof.dat');

%%make it interesting
hn1 = floor(nx/3*2); 
hn2 = floor(nx/5*4); 
hny = floor(ny/2);
ccof (hn1:hn2,hny) = 1;
%end interesting

%VARIABLES
vccof = ccof; %% volumetric version of ccof
for i=1:nx
    for j=1:ny
        for k=1:nz
            vccof(i,j,k) = ccof(i,j,k)*dx(i)*dy(j)*dz(k);
        end
    end
end

%DATA
%% initial condition set-up or read from file
nsolprev = zeros(nx,ny,1);
midpointx = floor([nx/3*2,nx/6*5]);
midpointy = floor([ny/3*2,ny/6*5]);
nsolprev(midpointx(1):midpointx(2) ,midpointy(1):midpointy(2) ,1) = 1;

%DISPLAY
if ifshow
    clf; 
    %plot(xc,nsolprev);
    contourf(xc,yc,nsolprev');
    colorbar;
end

%
%ALGORITHM
%% choose dt based on CFL condition constraint
datcfl = vccof./umag; 
dtcfl = min(datcfl(1:end));

if dt0>0
    dt = min(0.5*dtcfl,dt0);
    %dt = min(dtcfl,dt0);
else
    dt = 0.5*dtcfl;
    dt = dtcfl;
end

fprintf('dt cfl=%g dtrequested=%g dt Set to %g\n',dtcfl,dt0,dt);
pause;
%%
tsteps=0:dt:Tend; 
nsteps = length(tsteps); % approximate length of time steps
%%%
fprintf('Advection solver on domain [%g,%g] x [%g,%g] of dim %d x %d \n',x0,xend,y0,yend,nx,ny);
pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n =1; 
t = 0;
while 1
    tprev = t;
    t = t + dt; 
    if t> Tend
        t=Tend; 
        dt = t-tprev; 
    end
    tsteps(n) = t;
%
[concnew,total,inflow,outflow] = advect_step(tprev,t,vx,vy,vz,vccof,nsolprev,void,xc,yc,zc,dx,dy,dz);
nsol = concnew;
totalc = sum(total(1:end)); maxc = max(nsolprev(1:end));
%% post-process and display solution
if ifshow > 0 && (mod(n,ifshow) == 0 || n==1 || n==nt)
showtran(xc,yc,zc,xxc,yyc,nsol(:,:,1),void,t,n,totalc);
fprintf('t=%g nt=%d totalc=%g maxc=%g \n',t,n,totalc,maxc);
pause(0.05);
end

nsolprev = nsol;
%%
if t>= Tend
    break; 
else
    n = n +1; 
end
end %% time stepping loop
%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% actual advection time step
function [conc,totalconc,inflow,outflow] = advect_step(told,tnew,vx,vy,vz,phi,concold,void,xc,yc,zc,dx,dy,dz)
%%
%% take advection time step from <told> to <tnew> with fluxes in <vx,vy,vz>
%% previous time step value in <concold>, keyout <void>
%% Dirichlet value on the boundary if needed taken using dirichlet_fun, and <xc,xy,zc,dx,dy,dz>
%%
totalconc = concold.*phi; %% old time step values multipled by volumetric ccof given by <phi>
dt = tnew - told; 
t = tnew;

[nx,ny,nz] = size(concold);

if nargout>2
    inflow=0;
    outflow=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% advect  in x - direction
for j=1:ny
    for k=1:nz
        for l=1:nx+1 %% loop over the edges
            quant = 0;
             %% find the quantity to be advected,possibly from bdary
             if l>=2 && l<= nx
                 if vx(l,j,k) >= 0
                     quant = concold(l-1,j,k);
                 else 
                     quant = concold(l,j,k); 
                 end
             elseif l == 1
                 if vx(l,j,k) >= 0
                     dirval = dirichlet_fun (xc(l) - dx(l)/2, yc(j),zc(k), t);
                     quant = dirval;
                 else
                     quant = concold(l,j,k); 
                 end
             else % l == nx+1
                if vx(l,j,k) >= 0
                    quant = concold(l-1,j,k);
                else
                    dirval = dirichlet_fun (xc(l-1) + dx(l-1)/2, yc(j),zc(k),t) ;
                    quant = dirval;
                end
             end
            %% now that we know {quant}, update {totalconc}
            if l<=nx
                if void(l,j,k) == 1
                    totalconc(l,j,k) = totalconc(l,j,k) + dt*vx (l,j,k)*quant; 
                end
            end
            if l>1 %% or boundary condition
                if void(l-1,j,k)
                    totalconc(l-1,j,k) = totalconc(l-1,j,k) - dt*vx (l,j,k)*quant;
                end
            end

            %% if of interest, calculate influx or outflux across the boundary
            if l ==1 && void(l,j,k)==1
                inflow = inflow + dt*vx(l,j,k)*quant;
            end
            if l ==nx +1 && void(l-1,j,k)==1
                outflow = outflow - dt*vx(l,j,k)*quant;
            end

        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% advect in y - direction
for k=1:nz
    for j=1:nx
        for l=1:ny+1 %% loop over the edges
            quant=0;
            %% find the quantity to be advected,possibly from bdary
            if l>=2 && l<= ny
                if vy(j,l,k) >= 0
                    quant = concold(j,l-1,k);
                else
                    quant = concold(j,l,k); 
                end
            elseif l == 1
                if vy(j,l,k) >= 0
                    dirval = dirichlet_fun (xc(j),yc(l) - dy(l)/2,zc(k), t);
                    quant = dirval;
                else
                    quant = concold(j,l,k); 
                end
            else % l == nx+1
                if vy(j,l,k) >= 0
                    quant = concold(j,l-1,k);
                else
                    dirval = dirichlet_fun (xc(j), yc(l-1) + dy(l-1)/2,zc(k),t) ;
                    quant = dirval;
                end
            end
            %% now that we know {quant}, update {totalconc}
            if l<=ny
                if void(j,l,k) == 1
                    totalconc(j,l,k) = totalconc(j,l,k) + dt*vy (j,l,k)*quant; 
                end
            end
            if l>1 %% or boundary condition
                if void(j,l-1,k)
                    totalconc(j,l-1,k) = totalconc(j,l-1,k) - dt*vy (j,l,k)*quant;
                end
            end

            %% if of interest, calculate influx or outflux across the boundary
            if l ==1 && void(j,l,k)==1
                inflow = inflow + dt*vy(j,l,k)*quant;
            end
            if l ==ny +1 && void(j,l-1,k)==1
                outflow = outflow - dt*vy(j,l,k)*quant;
            end

        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% advect in z - direction
%% post-process: bring back the unknown.
conc = totalconc ./phi;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function showtran(xc,yc,zc,xp,yp,conc,void,tstep,nt,totsum)

%% display the numerical solution
vconc = conc(:,:,1); 
vconc(void==0) = nan;
%plot(xc,vconc);
contourf(xc,yc,vconc');
colorbar;
%caxis([0,1]);
title(sprintf('t=%6.4f min=%6.4f max=%6.4f tot=%6.4f\n',tstep,min(conc(1:end)),max(conc(1:end)),totsum));
xlim([-1 1]);
ylim([-1 1]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dx,x,xc,nx,h] = setup_grid1d(nxdx,x0,xend)

%% if corner points provided, ignore x0, xend
if length(nxdx)>1
    x = nxdx; nx=length(x)-1; 
    dx=diff(x);
else
    nx = nxdx; 
    if nargin<2
        x0=0;
        xend=1;
    end
    dx =(xend-x0) / nx*ones(1,nx);
    x = 0*dx; x(1)=x0; 
    for j=2:nx+1 
        x(j)=x(j-1)+dx(j-1);
    end
end

% setup the cell centers in xc
xc = x(1:nx)+dx/2; 
h = max(dx);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function value = dirichlet_fun(x, y, z, t)

if x>0.5 && y>0.5
    value = 0; 
else   
    value = 0;
end

end
