% SIR epidemiological model for MoHANET with 2-D Truncated Lévy Walk (TLW)
close all
clear 
clc

tic
%% Parameters
% SIR parameters
tau_SIR = 1; %time step for SIR model (day)
ts_SIR = 300; % total simulation time for SIR (day)
t = 0:tau_SIR:ts_SIR; %time vector for SIR model
gamma = 140; % threshold values
mu_R = [117 120 123]; % mean no of received droplets
sigma_R = 10; % std deviation of the no of received droplets

% TLW parameters
MC = 10^3; %no of Monte Carlo loops
no_of_nodes = 1000;
x_max = 2000; % x-dimension of the simulation area (m)
y_max = 2000; % y-dimension of the simulation area (m)
t_sim = 12*3600; % simulation time (s)
tau = 10; %time step (s)
d_inf = 1; % infection range (m)
no_of_steps = t_sim/tau; %no of steps
% theta = deg2rad(60); %beamwidth of the node
alpha_r = 1.6; % shape parameter of alpha-stable distribution for flight
c_r = 10; %scale parameter of alpha-stable distribution for flight
beta_r = 1; mu_r = 0; % skewness and location parameter - beta = 1, mu = 0 for Lévy distribution
alpha_t = 0.8; % shape parameter of alpha-stable distribution for time
c_t = 1; %scale parameter of alpha-stable distribution for time
beta_t = 1; mu_t = 0; % skewness and location parameter - beta = 1, mu = 0 for Lévy distribution
tau_p = 1000; %truncation limit for pause time distribution (s)
r_1 = 1000; %truncation limit for flight length distribution (m)

%% Determine the number of infected nodes with Truncated L�vy Walk (Rhee et al, 2011)
% parpool('local',7) %starts the parallel pool 
% parfor i_mc = 1:MC %parallel computing - If you activate parallel
% computing then comment the line below
for i_mc = 1:MC
    % Determine the positions of the nodes
    Nx = zeros(no_of_nodes,no_of_steps); % initialize x-position of the nodes (m)
    Ny = zeros(no_of_nodes,no_of_steps); % initialize y-position of the nodes (m)
    Nx(:,1) = x_max.*rand(1,no_of_nodes) ; % choose a random initial x-position with U(0,x_max)
    Ny(:,1) = y_max.*rand(1,no_of_nodes) ; % choose a random initial y-position with U(0,y_max)
    for i = 1:no_of_nodes
        phi = 2*pi.*rand(1,1*no_of_steps); %determine random direction with U(0,2*pi)
        delta_r_temp = random('Stable',alpha_r,beta_r,c_r,mu_r,1,1*no_of_steps); %random flight length from alpha-stable distribution
        delta_r = delta_r_temp(delta_r_temp > 0 & delta_r_temp <= r_1); %truncate
        delta_t_p_temp = random('Stable',alpha_t,beta_t,c_t,mu_t,1,1*no_of_steps); %random pause time from alpha-stable distribution
        delta_t_p = delta_t_p_temp(delta_t_p_temp > 0 & delta_t_p_temp <= tau_p); %truncate
        delta_r_temp = []; delta_t_p_temp = []; %clear variables from memory
        k = 1; k_r = 1;
        while (k <= no_of_steps) 
            % Flight time calculation
            if delta_r(k_r) < 500
                k_f = 30.55; rho_f = 0.89; %fitting parameters for flight time
            else
                k_f = 0.76; rho_f = 0.28; %fitting parameters for flight time
            end
            delta_t_f = k_f*delta_r(k_r)^(1-rho_f); %flight time  
            n_f = round(delta_t_f/tau); %calculate how many steps this flight takes      
            delta_x = delta_r(k_r)*cos(phi(k_r)); % Step length on x-axis 
            delta_y = delta_r(k_r)*sin(phi(k_r)); % Step length on y-axis       
            Nx(i,k+1:k+n_f) = Nx(i,k) + ((delta_x/n_f):(delta_x/n_f):delta_x); % update the position of the i^th node on x-axis for the flight
            Ny(i,k+1:k+n_f) = Ny(i,k) + ((delta_y/n_f):(delta_y/n_f):delta_y); % update the position of the i^th node on y-axis for the flight

            % Pause time
            n_p = round(delta_t_p(k_r)/tau); %calculate how many steps this flight takes            
            Nx(i,k+1+n_f:k+n_f+n_p) = Nx(i,k+n_f); % update the position of the i^th node on x-axis along pause time
            Ny(i,k+1+n_f:k+n_f+n_p) = Ny(i,k+n_f); % update the position of the i^th node on y-axis along pause time

            % Border control of the nodes: If the next step of the node is
            % outside of the simulation zone, then the node is boncuking from the border.
            if max(Nx(i,k+1:k+n_f+n_p)) > x_max  
                i_temp = find(Nx(i,:) > x_max, 1); 
                Nx(i,i_temp:k+n_f+n_p) = x_max - mod(Nx(i,i_temp:k+n_f+n_p), x_max); % bouncing
            elseif min(Nx(i,k+1:k+n_f+n_p)) < 0
                i_temp = find(Nx(i,:) < 0, 1); 
                Nx(i,i_temp:k+n_f+n_p) = x_max - mod(Nx(i,i_temp:k+n_f+n_p), x_max); % bouncing                
            end

            if max(Ny(i,k+1:k+n_f+n_p)) > y_max  
                i_temp2 = find(Ny(i,:) > y_max, 1); 
                Ny(i,i_temp2:k+n_f+n_p) = y_max - mod(Ny(i,i_temp2:k+n_f+n_p), y_max); % bouncing
            elseif min(Ny(i,k+1:k+n_f+n_p)) < 0
                i_temp2 = find(Ny(i,:) < 0, 1); 
                Ny(i,i_temp2:k+n_f+n_p) = y_max - mod(Ny(i,i_temp2:k+n_f+n_p), y_max); % bouncing                
            end
            k = k+n_f+n_p;
            k_r = k_r+1;    

%                 Nx(i,k+1:k+n_f+n_p) = 0; Ny(i,k+1:k+n_f+n_p) = 0; %remove the values out of borders
%                 k_r = k_r+1;
%                 continue; %exit loop and repeat the same step
%             k = k+step_temp;
        end
        delta_r = []; delta_t_p = [];
    end
    Nx = Nx(:,1:no_of_steps); Ny = Ny(:,1:no_of_steps); %truncate the residues
%     step_len_tot(i_mc,:) = reshape(step_len,[1,size(step_len,1)*size(step_len,2)]);%total step lengths for each MC loop
    %% Determine the distances and infection state of the nodes
    N_inf = zeros(no_of_nodes, no_of_steps); %infection state of the nodes at each time step
    N_inf(1,:) = 1; %Node 1 is initially infected

    N_inf = infection2(no_of_nodes, no_of_steps, d_inf, Nx, Ny, N_inf);% Trace the infection state of the nodes - omnidirectional
%     N_inf = infection_bw(no_of_nodes, no_of_steps, d_inf, Nx(:,1:no_of_steps), Ny(:,1:no_of_steps), N_inf, theta, phi(:,1:no_of_steps));% Trace the infection state of the nodes with beamwidth - directional
    no_inf_nodes(i_mc,:) = sum(N_inf);% Total number of infected nodes at each time step
%     Nx = []; Ny=[];  N_inf = [];%clear variables from memory
end

N_c = mean(no_inf_nodes); % average no of infected nodes at each time step

N_c_avg = mean(N_c)/(t_sim/(24*3600)); %Average contact rate (per day)

%% Determine the infection rate
P_inf = qfunc( (gamma - mu_R) / sigma_R );

for i = 1:length(P_inf)
   beta(i) = P_inf(i) * N_c_avg; %average infection rate
end


%% SIR Model
I(1:length(mu_R),1) = 1; % no of infected humans
S(1:length(mu_R),1) = no_of_nodes - I(1:length(mu_R),1); % number of susceptible humans
R(1:length(mu_R),1) = 0; % no of recovered humans
beta_2 = 0.025; %Recovery rate

% With average infection rate
for i = 1:length(mu_R)
    for j = 1:length(t)-1
        dS = (- beta(i)*I(i,j)*S(i,j))*tau_SIR/no_of_nodes;
        S(i,j+1) = S(i,j) + dS;
        dI = (beta(i)*I(i,j)*S(i,j)/no_of_nodes - beta_2*I(i,j))*tau_SIR;
        I(i,j+1) = I(i,j) + dI;
        dR = (beta_2*I(i,j))*tau_SIR;
        R(i,j+1) = R(i,j) + dR;
    end
end

h = figure;
plot(t, I(1,:),'b--', t, I(2,:),'r:', t, I(3,:),'g-', 'LineWidth',1.25); 
legend(['\mu_R = ', num2str(mu_R(1))], ['\mu_R = ', num2str(mu_R(2))], ['\mu_R = ', num2str(mu_R(3))]);
xlabel('Time (days)'); ylabel('Total number of actively infected humans');
grid on;
ylim([0 no_of_nodes]);

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
% print(h,sprintf('SIR_TLW_mu_R_plot_MC_%d.pdf', MC),'-dpdf','-r0')%save as pdf 
% savefig(h,sprintf('SIR_TLW_SIR_plot_dinf_%.1f_alpha_r_%.1f_MC_%d.fig',d_inf,alpha_r,MC)); %save the figure file
toc


