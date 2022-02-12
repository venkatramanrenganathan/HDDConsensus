%% NECSYS 2022 Simulation Code
% 
% This code simulates the History data-driven consensus protocol. The code
% is tested with Matlab 2020b on a MacOS. 
% 
% Date: 10 February, 2022 
%
% Contributing Authors:
% 1) Angela Fontan - Email: angela.fontan@kth.se
% 2) Venkatraman Renganathan - Email: venkat@control.lth.se
%
% (C) Angela Fontan, Venkatraman Renganathan, 2022. 
%
% This program is a free software: you can redistribute it and/or modify it
% under the terms of the GNU lesser General Public License, either version 
% 2018b, or any later version. This program is distributed in the hope that 
% it will be useful, but WITHOUT ANY WARRANTY. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Set some default plot specifications
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaulttextInterpreter','latex');

% clear the variables and command history
close all; clear variables; clc;

% Flag to decide to load a presaved data or not
loadData = 0; 

% Flag to decide to save the data
saveFig = 1;

% Flag to decide loading a predefined initial conditions (history) 
loadIC = 1;

% Create the history and perform HDD consensus protocol
if ~loadData
    
    % Define the protocol for non-cooperative agents
    f_malicious = @(idx_m,it) [2.5+0.3*randn;-3*0.5*randn;-2*0.1*randn];

    % Define the connection graph of agents
    nr = 10;            % Number of cooperative agents    
    nm = 3;             % Number of non-cooperative agents 
    n = nr+nm;          % Total number of agents
    idx_r = 1:nr;       % Indices of cooperative agents
    idx_m = nr+(1:nm);  % Indices of non-cooperative agents
    pr = 0.4;           % Probability of an edge between cooperative agents
    fNi = 0;
    while ~fNi
        Ni = repmat(1:n,n,1);
        WNi = triu(rand(n).*((pr-rand(n)) > 0)) + diag(rand(n,1)); 
        WNi(:,idx_m) = 1;
        WNi = WNi + WNi.';
        Ni = Ni.*(WNi > 0);
        fNi = 1;
        %fNi = nnz(conncomp(graph(WNi))>1)==0;
    end
    
    % Check & display if the resulting graph is connected or not.
    fprintf('Graph of cooperative agents is connected? %d\n',nnz(conncomp(graph(WNi(idx_r,idx_r)))>1)==0);    

    % create graph and history
    t0 = -15;   % Starting point(index) of the history
    tf = 200;   % End point of Simulation 
    t = t0:tf;  % Time index vector
    Neighbors = repmat({Ni}, 1, tf - t0 + 1); % Same set of neighbors at all time steps
    
    % If loadIC = True, preload from mat file, else create the history
    if ~loadIC
        x_history = 2.2*randn(n,0-t0+1)+0.1*(1:n)';
        save('ic.mat','x_history');
    else
        load('ic.mat');
    end

    % --- History-based consensus
    par_HBC.T = 15;   % History length 
    emax = 0.5;       % maximum (epsilon) confidence bound
    emin = 0.01;      % minimum (epsilon) confidence bound
    
    % Generate a decreasing sequence of random epsilons from U[emin, emax]
    par_HBC.epsilon = sort(emin+(emax-emin)*rand(1,length(t)), 'descend');
    
    % Vector of discount factor
    nu = 0.05:0.05:0.95;
    
%     % Following three parameters is for Antoine Girard's paper
%     eps_R = 1; eps_rho = .9; par_HBC.epsilon = eps_R*eps_rho.^t;   
    
    % Define a struct to hold simulation variables
    comparison.trust = cell(1,length(nu)); 
    comparison.x = cell(1,length(nu)); 
    comparison.W = cell(1,length(nu));
    comparison.x_end = zeros(nr,length(nu));

    % Loop across different discount factor values
    for inu = 1:length(nu)
        
        % Cell to hold the weights for all time steps for inu
        W_t = cell(1,tf-t0+1);
        
        % Set the discount factor as inu
        par_HBC.nu = nu(inu);
        
        % Call & execute the history based consensus for this time step
        [x_HBC,W_t_HBC,trust] = f_HBC(x_history,t,idx_r,idx_m,Neighbors,W_t,f_malicious,par_HBC);
        
        % Store the results in respective cells
        comparison.W{inu} = W_t_HBC;
        comparison.x{inu} = x_HBC;
        comparison.trust{inu} = trust;
        comparison.x_end(:,inu) = x_HBC(idx_r,end);        
    end
    
    % Save all the data and close any existing figures
    save('Data_ex_nu.mat'); close all;          
    
    % Call the plotter function
    f_fig_nu(saveFig,t,idx_r,idx_m,nu,comparison,par_HBC);
    
else
    % Just plot the simulated data 
    % load the simulated data close any existing figures
    load('Data_ex_nu.mat'); close all;     
    
    % Call the plotter function
    f_fig_nu(saveFig,t,idx_r,idx_m,nu,comparison,par_HBC);
    
end