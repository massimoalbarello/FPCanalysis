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

% Open-Loop System
average_A_P = average_A - B*Kp*C;

% Definition of variables for storage of state evolution and Measurements
% evolution
x_k_average_P = zeros(n , t_end+1); 
y_k_average_P = zeros(n_selfish , t_end+1); 

x_k_average_P(: , 1) = x_0;
y_k_average_P(: , 1) = C * x_0;

% Building the sequence of state matrices for P-Controlled Closed Loop 
for k  = 1:t_end
    %A_cl_sequence(: , : , k) = [A_sequence(: , : , k)-B*Kp*C , B*Ki ; -C , 1];
    %x_k_PI(: , k+1) = sat_function(A_cl_sequence(: , : , k) * x_k_PI(: , k) + [B*Kp]*ref);
    x_k_average_P(: , k+1) = sat_function(average_A_P * x_k_average_P(: , k) + [B*Kp]*ref);
    y_k_average_P(: , k+1) = C * x_k_average_P(: , k+1);
end

figure(101) ;  hold on;
plot(0:1:t_end , x_k_average_P(1:3 ,:) ,  'LineWidth' , 1.5); hold on;
plot(0:1:t_end, x_k_average_P(n_selfish+1 ,:),  'LineWidth' , 1.5); 
plot(0:1:t_end , y_k_average_P(1 ,:) , 'LineWidth' , 1.5);
plot(0:1:t_end, ref_seq, 'k .' , 'MarkerSize' , 1.1);
legend( 'Coordinator 1' ,'Coordinator 2' ,'Coordinator 3' , 'Standard Agent 1' , 'Network average' , 'Reference');
% title('Global, saturation, Mean reference, P');
hold off;