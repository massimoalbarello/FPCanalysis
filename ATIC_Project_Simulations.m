%% RANDOMIZED MODEL ON A COMPLETE GRAPH
% Defining the Adjacency matrix dimensions (n nodes)
n = 50;
% Defining the number of nodes to be sampled
p = 5;
%Fix simulation length
t_end = 500;       
%plotting variable
plotting = false;

% Defining the initial conditions on the opinions
x_0 = zeros(n , 1);
% for i=1:n
%     x_0(i , 1) = round(rand(1));
% end

%Seed for initial conditions
x_0 = [1;0;1;0;1;0;1;0;0;1;0;0;1;1;1;1;1;1;0;1;0;0;1;1;1;1;0;1;1;1;0;1;0;0;1;0;0;1;1;0;1;0;1;0;1;1;0;1;0;1];
    

%Computing the random masks sequence and A(k)
mask_sequence = zeros(n , n , t_end);
for i = 1:t_end
    mask_sequence(: , : , i) = mask_samples(n , p);
end
A_sequence = mask_sequence * (1/(p+1));


% Computing the average of the sequence of masks to check if correct
% distribution
average_mask = p/(n-1)*(ones(n)-eye(n))+eye(n);
average_A = average_mask * (1/(p+1));

%Check row stochasticity of A (satisfied)
row_sums = sum(average_A,2);
col_sums = sum(average_A,1);
row_sum = sum(A_sequence,2);
col_sum = sum(A_sequence,1);

%State evolution of the dicrete time system (Randomized)
x_k = zeros(n , t_end); 
x_k(: , 1) = x_0;           %Initial state = Initial conditions
for k = 1:t_end
    x_k(: , k+1) = A_sequence(: , : , k) * x_k(: , k);
end
% x_k;

%State evolution of the dicrete time system (Expected Value of the Sequence)
x_k_expected = zeros(n , t_end+1);
x_k_expected(: , 1) = x_0;           %Initial state = Initial conditions
for k = 1:t_end
    x_k_expected(: , k+1) = average_A * x_k_expected(: , k);
end
% x_k_expected;

%Plotting the DTime evolution of all the states (Comparing with average
%matrix)

if(plotting == true)
    figure(1); plot(0:1:t_end , x_k , 'LineWidth' , 1.1); grid on;
    figure(2); plot(0:1:t_end , x_k_expected, 'LineWidth' , 2); grid on;
end

% Computing the Laplacian matrix associated to the sequence of A(k)
L_sequence = zeros(n , n , t_end);
for i = 1:t_end
    L_sequence(: , : , i) = eye(n)- A_sequence(: , : , i);
end
average_L = eye(n) - average_A;

%Plotting the evolution using the sequence of laplacian matrices
step = 0.01;
counter = 1;
t_span = 0:step:t_end;
x_t = zeros(n , length(t_span));
x_t(: , 1) = x_0;
x_init = x_0;
for i = 1:t_end
    for k = 1:(1/step)
        x_t(: , (1/step)*(i - 1)+ k + 1) = expm(-L_sequence(: , : , i)*t_span(k))*x_init;
    end
    x_init = x_t(: , 100*(i - 1)+ k + 1);
end

if(plotting == true)
    figure(3); plot(t_span , x_t( : ,:), 'LineWidth' , 1.1); hold on; grid on;
    % plot(0:1:t_end , x_k(1:2 ,:), 'LineWidth' , 1.3)
end

%Computing the solution for the continuous time system (Expected L)
x_t_expected = zeros(n , length(t_span));
for k = 1:length(t_span)
    x_t_expected(: , k) = expm(-average_L*t_span(k))*x_0;
end


%Plotting the state evolution in continuous time (for expected L)
if(plotting == true)
    figure(4); plot(t_span , x_t_expected, 'LineWidth' , 2); grid on;
end



%% IMPLEMENTING THE FEEDBACK CONTROL LAW

%% STEP 1
% Assumptions:
% - One selfish agent IOTA NODE
% - Randomized Adjacency
% - Myopic Selfish Agent C not strictly positive C = [c_1 , ... , c_n] only
%   m c_i's are non-zero, with m<n
% - complete undelying network
% - The mean has to converge to a reference ref

n_selfish = 1;

%Define reference value
ref=mean(x_0);

% Define the closed loop system state space representation (A , B , C , D)
B = [ones(n_selfish , 1) ; zeros(n - n_selfish , 1)];
C = (1/n)*[ones(1 , n)];
Kp = 10;
Ki = 0.01;

%Initial condition for the augmented states
x_0=[x_0;0];
x_k_cloop = zeros(n + n_selfish , t_end+1);
x_k_cloop(: , 1) = x_0;

A_cl_sequence = zeros(n + n_selfish , n + n_selfish , t_end);
for k  = 1:t_end
    A_cl_sequence(: , : , k) = [A_sequence(: , : , k)-B*Kp*C , B*Ki ; -C , 1];
    x_k_cloop(: , k+1) = A_cl_sequence(: , : , k) * x_k_cloop(: , k) + [B*Kp ; 1]*ref;
end
figure(5) ; plot(0:1:t_end , mean(x_k_cloop,1) ,  'LineWidth' , 2); grid on; hold on; legend();
plot(0:1:t_end , x_k_cloop(51 , :) ,  'LineWidth' , 2); hold off;


%% STEP 2
% Assumptions:
% - One selfish agent IOTA NODE
% - Randomized Adjacency
% - Myopic Selfish Agent C not strictly positive C = [c_1 , ... , c_n] only
%   m c_i's are non-zero, with m<n
% - complete undelying network
% - The error wrt the other nodes has to converge to a reference ref = 0

% Define the closed loop system state space representation (A , B , C , D)
B = [ones(n_selfish , 1) ; zeros(n - n_selfish , 1)];
C = [n-n_selfish , - ones(1 , n-n_selfish)];
Kp = 0.01;
Ki = 0.001;

% Defgining new reference
ref = 0.001;

A_cl_sequence = zeros(n + n_selfish , n + n_selfish , t_end);
for k  = 1:t_end
    A_cl_sequence(: , : , k) = [A_sequence(: , : , k)-B*Kp*C , B*Ki ; -C , 1];
    x_k_cloop(: , k+1) = A_cl_sequence(: , : , k) * x_k_cloop(: , k) + [B*Kp ; 1]*ref;
end
figure(6) ; plot(0:1:t_end , x_k_cloop(:,:) ,  'LineWidth' , 2); grid on; hold on; legend();
plot(0:1:t_end , x_k_cloop(51 , :) ,  'LineWidth' , 2); hold off;
