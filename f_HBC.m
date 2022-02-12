function [x,W_t,trust] = f_HBC(x_history,t,idx_r,idx_m,Neighbors,W_t,f_malicious,par_HBC)

% Unbox the parameters
nu = par_HBC.nu;                   % Discount factor
T  = par_HBC.T;                    % Length of History
epsilon = par_HBC.epsilon;         % Confidence bound sequence
n = length(idx_r) + length(idx_m); % Total number of agents
trust.delta = cell(1,length(t));   % Store the discounted importance vector 
trust.mu = cell(1,length(t));      % Store the estimated mean

% Define placeholder for states
x = [x_history, zeros(n,length(t) - size(x_history,2))]; 

% Loop through till the final time
for it = find(t==0):(length(t)-1)
    
    % Time index sequence $\kappa_{t,T}$ starting from any "t"
    k_t_T = (t(it) - T + 1):(t(it)); 
    
    % Indices denoting $\kappa_{t,T}$
    idx_k_t_T = (it - T + 1):it;   
    
    % Extract the epsilon sequence
    e = epsilon(idx_k_t_T);    
    
    % Set membership Estimation parameters
    W_t{it}   = zeros(n);    % Weights at time "it"
    mu_i      = cell(1,n);   % estimated mean at time "it"
    N_i_k     = cell(1,n);   % Set membership counter at time "it"
    delta_i_j = cell(1,n);   % Discounted importance vector
    
    % Loop through all cooperative agents 
    for i_agent = idx_r
        
        % Find neighbors whose values are less than epsilon_k
        states_less_e = abs(x(setdiff(1:n,i_agent),idx_k_t_T)-x(i_agent,idx_k_t_T))-repmat(e,n-1,1)<0;
        
        % Define variables 
        count = 0;  
        N_i   = zeros(n,T);
        
        % Count the number of people in the ball with radius epsilon_k
        for k = idx_k_t_T
            count = count + 1; 
            N_i(:,count) = Neighbors{k}(i_agent,:)';
        end
        N_i_0 = N_i; 
        N_i_0(i_agent,:)=[]; 
        N_i_k{i_agent} = N_i_0.*states_less_e;
        
        % Compute the discounted importance vector 
        delta_i_j{i_agent}=(N_i_k{i_agent}>0)*diag(nu.^(t(it)-k_t_T));
        
        % Define the estimated mean placeholder for ith cooperative agent
        mu_i{i_agent} = zeros(n,1); 
        mu_i{i_agent}(i_agent) = 1;
        
        % Compute the estimated mean
        mu_i{i_agent}(setdiff(1:n,i_agent)) = sum(delta_i_j{i_agent},2)/T; % this is zeta basically
        
        % Compute the normalized weight
        W_t{it}(i_agent,:)=mu_i{i_agent}'/norm(mu_i{i_agent},1);
    end
    
    % Prepare the output to be returned
    trust.mu{it}    = mu_i;
    trust.delta{it} = delta_i_j;
    
    % Define the non-cooperative agents values using its protocol "f_malicious"
    x(idx_m,it+1) = f_malicious(idx_m,it);
    
    % Execute the HDD Consensus protocol with the above weights
    x(idx_r,it+1) = W_t{it}(idx_r,idx_r)*x(idx_r,it) + W_t{it}(idx_r,idx_m)*x(idx_m,it);
    
end


end

