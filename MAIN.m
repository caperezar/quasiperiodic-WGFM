% Written by Carlos Perez-Arancibia (caperezar@gmail.com)

close all
clear

%% Physical parameters
angle = -pi/2;    % angle of incidence with respect to the horizontal
L     = 2;        % array period
pol   = 'TM';     % polarization

k1    = 12;        % exterior wavenumber
k2    = 16;        % interior wavenumber

alpha =  k1*cos(angle);
beta  = -k1*sin(angle);

lambda = 2*pi\k1;

% the physical date is stored in this data structure:
pd    = physical_data([k1 k2]',pol,angle,L);


%% Parameters of the discretization:
c     = 0.5;       % parameter that controls the decay of the window function
tol   = 1e-2;      % parameter to control when to use L'Hopital rule

delta = 3/4*k1;    % parameter that selects the correction modes 
A     = 5*lambda; % window size

h_corr   = c*A;         % height used for the correction term
h_energy = 0.9*h_corr;  % height used for the energy test (error)

np_obs   = 200;    % number of discretization points on the obstacle's boundary 
np_walls = 400;   % number of discretization points on the cell walls
np_hrz   = 10;   % number of points used in the discretization of the horizontal lines      

%% Geometry:
geo = geometry; 

num_obstacles  = 1;             % number of obstacles
obstacle_type  = 'kite';        % obstacle type

num_cell_walls = 2;             % number of cell walls
cell_wall_type = 'line_segment';  % curve type 

geo = geo.cctor(num_obstacles,obstacle_type,num_cell_walls,cell_wall_type);

% Construction of the obstacle(s) discretization (\Gamma_1):
rad_obs  = 1/2;    % radius of the obstacle
cntr_obs = [0 0];  % centre of the obstacle
geo.obstacles = geo.obstacles.cctor(floor(np_obs/2),cntr_obs,rad_obs);

% Construction of the discretization of the cell walls (\Gamma_2 and \Gamma_3):
start_left  = [-L/2 -A]; % starting point of left wall
end_left    = [-L/2 A];  % ending point of left wall
geo.segments(1)  = geo.segments(1).cctor(floor(np_walls/2),start_left,end_left,c);

start_right = [ L/2 -A]; % starting point of right wall
end_right   = [ L/2 A];  % ending point of right wall
geo.segments(2)  = geo.segments(2).cctor(floor(np_walls/2),start_right,end_right,c);

% Plot of the geometry:
geo.plot(1); hold off

%% Linear system matrices:

% main system matrices:
[M,W_A,E] = create_matrices_periodic(pd,geo);

% correction matrices:
[M_corr,Lp,Lm,dLp,dLm] = correction_matrix(pd,geo,h_corr,np_hrz);

% right hand side of the system
f_inc  = right_hand_side(pd,geo);

% linear system solution:
sol = (E+M*W_A+M_corr)\f_inc;
% [sol,FLAG,RELRES,ITER] = gmres(E+M*W_A+M_corrr,f_inc,1000,1e-6,1000,E);

% turning linear system solution into densities:
dens = densities;
dens = dens.cctor(sol,pd,geo);

%% Conservarion of energy error:
[Error,R,T] = energy_test(h_energy,geo,pd,dens,Lp,Lm,dLp,dLm);
% R: reflection coefficient:
% T: transmission coefficient
Error

%% Plotting the solution:


h = 0.05;         % grid size
x = -L/2:h:L/2;   % x-grid
y = -1.5*h_energy:h:1.5*h_energy; % y-grid

option = 't'; % generate the total field
n_plot = 2;   % plot the field on the cless -n_plot:n_plot


[U,X,Y] = evaluate_field_super_cell(x,y,pd,geo,dens,Lp,Lm,dLp,dLm,option);

% plot total field:
figure(2);hold on

for c=-n_plot:n_plot
    surf(X+c*L,Y,real(U*pd.gm^c),'EdgeColor','none');
    for o=1:geo.n_o
        pts= [geo.obstacles(o).x(:,1) geo.obstacles(o).x(:,2) 100*ones(size(geo.obstacles(o).x(:,1)))];
        plot3(pts(:,1)+c*L,pts(:,2),pts(:,3),'k','LineWidth',1)
    end
end

view(2)
caxis([-2 2])
colormap(brewermap([],'*RdBu'));
axis equal
axis tight
axis off



