% Defining the Adjacency matrix dimensions (n nodes)
n = 100;
% Defining the number of nodes to be sampled
p = 3;
%Fix simulation length
t_end = 50; 

%Seed for initial conditions in which we have a clear majority for op. 1
x_0 = [0;1;1;1;1;1;1;1;0;1;1;1;1;0;0;1;1;1;0;1;1;1;0;1;1;1;0;1;1;1;0;1;1;1;0;0;0;0;1;1;0;0;1;0;1;0;1;1;1;1;1;1;1;1;1;1;1;0;1;1;1;0;1;0;1;0;0;1;1;1;0;1;1;1;0;1;1;0;1;1;1;0;1;1;0;1;1;1;1;1;0;1;0;1;1;0;1;1;1;1];

% Reference and Coordinating Nodes initialization
% number of coordinators
n_selfish = 5;

% Refererence chosen to be the mean of the standard agents
ref=mean(x_0(n_selfish+1:end))*ones(n_selfish , 1);

% Reference sequence (for plots)
ref_seq = mean(x_0(n_selfish+1:end)) * ones(t_end+1 , 1);

% Computing the average Mask
average_mask = p/(n-1)*(ones(n)-eye(n))+eye(n);

% State-Space Representation 
average_A = (1/(p+1))*average_mask;
B = [eye(n_selfish) ; zeros(n - n_selfish , n_selfish)];
C = (1/n)*[ones(n_selfish , n)];
Kp = 2*eye(n_selfish);
Ki = 0.01*eye(n_selfish);

% Definition of variables for storage of the augmented state evolution and the Measurements
% evolution
x_0_aug=[x_0; zeros(n_selfish,1)];
x_k_average_PI = zeros(n + n_selfish , t_end+1);
x_k_average_PI_sat = zeros(n + n_selfish , t_end+1);
y_k_average_PI = zeros(n_selfish , t_end+1);
x_k_average_PI(: , 1) = x_0_aug;

y_k_average_PI(: , 1) = C*x_0;

% Open-Loop System
average_A_PI = [average_A-B*Kp*C , B*Ki ; -C , eye(n_selfish)];
perturbation = [zeros(n_selfish , n) , Ki ; zeros(n , n+n_selfish)];
average_A_PI_np = average_A_PI-perturbation;

% Building the sequernce of state matrices for PI-Controlled Closed Loop 
for k  = 1:t_end
    x_k_average_PI(: , k+1) = average_A_PI * x_k_average_PI(: , k) + [B*Kp ; eye(n_selfish)]*ref;
    x_k_average_PI(1:n , k+1) = sat_function(x_k_average_PI(1:n , k+1));
    x_k_average_PI(n+1:end , k+1) = x_k_average_PI(n+1:end , k) + (ref - y_k_average_PI(: , k));
    y_k_average_PI(: , k+1) = C * x_k_average_PI(1:n , k+1);
end

%Plotting opinion Dynamics and integrated error
figure(103) ;  hold on;
plot(0:1:t_end , x_k_average_PI(1:3 ,:) ,  'LineWidth' , 1.5); hold on;
plot(0:1:t_end, x_k_average_PI(n_selfish+1:n_selfish+3 ,:),  'LineWidth' , 1.5);
plot(0:1:t_end, x_k_average_PI(n+1:n+3 ,:),  'LineWidth' , 1.5);
plot(0:1:t_end , y_k_average(1 ,:) , 'LineWidth' , 1.5);
plot(0:1:t_end, ref_seq, 'k .' , 'MarkerSize' , 0.7);
legend( 'Coordinator 1' ,'Coordinator 2' ,'Coordinator 3' , 'Standard Agent 1' , 'Standard Agent 2' , 'Standard Agent 3' , 'Integrated Error Coordinator 1' , 'Integrated Error Coordinator 2' ,  'Integrated Error Coordinator 3' , 'Network average' , 'Reference');
title('Global, saturation, Mean reference, PI');
hold off;