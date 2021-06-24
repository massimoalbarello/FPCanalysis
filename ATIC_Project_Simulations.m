%% RANDOMIZED MODEL ON A COMPLETE GRAPH
% Defining the Adjacency matrix dimensions (n nodes)
n = 100;
% Defining the number of nodes to be sampled
p = 5;
%Fix simulation length
t_end = 25;       
%plotting variable
plotting = false;
plotting_eig = false;
lap = false;
complete = false;
not_complete = true;

% Defining the initial conditions on the opinions
x_0 = zeros(n , 1);

% Random generator of initial condition using flip coin
for i = 1:n
    x_0( i , 1 ) = rand(1);
    if x_0(i , 1) > 0.25
        x_0(i , 1) = 1;
    else
        x_0(i , 1) = 0;
    end
end

%Seed for initial conditions in which we have a clear majority for op. 1
x_0 = [0;1;1;1;1;1;1;1;0;1;1;1;1;0;0;1;1;1;0;1;1;1;0;1;1;1;0;1;1;1;0;1;1;1;0;0;0;0;1;1;0;0;1;0;1;0;1;1;1;1;1;1;1;1;1;1;1;0;1;1;1;0;1;0;1;0;0;1;1;1;0;1;1;1;0;1;1;0;1;1;1;0;1;1;0;1;1;1;1;1;0;1;0;1;1;0;1;1;1;1];

%% GENERATING A COMPLETE TOPOLOGY AND SAMPLE FORM THE NEIGHBOURS, THE NODES ALWAYS SAMPLE THEMSELVES
if complete
    %Computing the random masks sequence and A(k)
    
    topology = ones(n);
    
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


%% GENERATING AN ARBITRARY TOPOLOGY (LOOK AT GRAPH_GEN FOR MORE DETAIL), THE NODES CAN SAMPLE FORM THE NEIGHBOURS ONLY, BUT ALWATS SAMPLE THEMSELVES
elseif not_complete
    
% Variable to check if the requirements on the topology are met
bad_topology = true;

while bad_topology
    % Generating a suitable underlying graph
    [topology] = gen_graph(n,p);

    %Checking irreducibility of the topology binary matrix
    power_A = eye(n);
    for i = 1:(n-1)
        power_A = power_A + topology^i;
    end
    
    Irreducibility = 1; %if one the matrix is irreducible, if zero it is reducible
    for i = 1:n
        for j = 1:n
            if power_A(i,j)== 0
                Irreducibility = 0;
            end
        end
    end

    % Checking the number of out neighbours including i itself
    sum_row = sum(topology , 2);
    Summation = 1;
    for i = 1:n
        if(sum_row(i)<(p+1))
            Summation = 0;
        end
    end

    
    if ((Summation==0) & (Irreducibility == 1))
        fprintf('Need more neighbours, Irreducible\n');
    elseif ((Summation==1) & (Irreducibility == 1))
        fprintf('Enough neighbours, Irreducible\n');
        topology = topology;
        bad_topology = false;
    elseif ((Summation==0) & (Irreducibility == 0))
        fprintf('Need more neighbours, Not Irreducible\n');
    elseif ((Summation==1) & (Irreducibility == 0))
        fprintf('Enough neighbours, Not Irreducible\n');
    else
         fprintf('I forgot one case\n');
         bad_topology = false;
    end
    
end

% Choose a network for the simulations from now on (play with it...)
topology = [1,1,1,1,1,1,0,0,1,0,1,1,1,1,0,1,1,1,1,1,0,0,0,0,1,0,1,0,0,1,1,1,1,0,0,1,0,1,1,1,1,1,0,1,0,0,0,1,1,0,1,0,1,1,0,1,0,1,0,0,1,1,1,1,1,0,1,0,0,0,0,0,1,1,1,0,0,1,0,1,1,0,0,0,1,0,0,1,0,0,0,1,1,0,1,0,0,1,0,0;1,1,1,1,0,0,1,0,1,0,1,0,0,1,0,0,0,1,0,0,0,0,0,0,1,1,0,1,1,0,1,0,1,1,0,0,0,1,1,0,1,0,1,1,0,0,0,1,0,1,0,0,1,0,0,1,0,1,0,0,1,0,0,0,0,0,0,1,1,0,1,0,0,1,1,0,0,1,1,0,0,0,1,0,1,1,1,0,1,1,1,0,0,1,1,1,1,1,0,0;1,1,1,0,1,1,0,1,1,1,1,1,0,0,0,1,1,0,1,1,1,0,1,0,1,1,1,1,0,1,0,1,0,0,0,0,1,0,1,1,1,1,1,1,0,0,1,1,0,1,1,1,0,0,0,1,0,0,1,1,0,0,1,1,0,0,0,0,0,1,1,0,0,0,0,0,1,1,1,1,0,0,0,1,0,1,1,1,0,0,1,0,0,1,1,0,0,1,0,0;1,1,0,1,1,0,0,1,1,0,1,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,1,1,1,0,1,1,0,0,0,0,0,1,1,0,0,1,1,0,0,0,0,1,0,0,0,1,1,0,0,0,1,1,0,1,1,1,0,1,0,0,1,1,0,0,1,0,0,0,1,0,1,0,0,0,1,1,1,0,1,0,1,0,1,0,1,1,0,0,0;1,0,1,1,1,0,1,1,1,1,1,0,0,1,0,0,1,1,1,0,1,0,0,1,1,1,0,1,0,0,0,1,0,0,0,1,1,0,0,0,1,0,1,1,1,1,0,0,1,0,1,0,0,1,0,0,1,1,0,0,0,0,1,0,0,1,0,1,1,1,0,0,1,1,0,1,1,1,1,0,0,1,0,0,0,0,0,1,0,0,1,1,0,0,1,1,1,1,1,0;1,0,1,0,0,1,0,1,1,0,1,1,1,0,1,1,0,1,0,1,0,1,1,0,0,0,1,0,0,1,0,1,1,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,1,1,1,0,1,0,1,1,0,0,0,0,0,1,1,0,0,1,1,0,1,1,0,0,0,0,0,1,0,0,0,1,1,1,0,1,1,0,1,1,0,0,0,1,0,0,1,1,0,0,0,0;0,1,0,0,1,0,1,0,0,1,1,1,0,1,0,1,0,1,0,0,0,1,0,0,1,0,0,1,1,0,1,1,0,1,1,1,0,0,0,0,0,0,0,1,1,1,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,1,1,0,1,0,0,0,0,0,1,0,0,0,1,0,0,0,1,1,0,0,0,0;0,0,1,1,1,1,0,1,1,1,1,0,0,1,1,1,0,1,0,0,1,1,1,0,1,1,0,1,0,0,0,0,0,0,0,1,1,0,0,1,1,0,0,0,1,1,1,1,0,0,0,0,1,0,0,0,0,1,1,0,0,1,0,0,1,1,0,0,0,1,0,0,0,1,0,1,1,1,1,0,0,1,1,0,1,0,0,0,0,1,1,0,1,1,1,0,0,0,0,0;1,1,1,1,1,1,0,1,1,1,0,0,1,0,1,0,0,1,0,0,0,1,1,1,0,0,0,0,0,1,0,0,1,0,1,1,1,1,0,1,0,1,0,0,1,0,0,0,0,1,1,0,1,0,1,0,0,1,1,0,0,0,0,1,1,1,0,1,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,0,0,1,1,0,0,0,0,1,1,1,1,1,0,0,0,1;0,0,1,0,1,0,1,1,1,1,0,1,0,1,1,1,0,1,0,1,1,0,1,1,0,0,1,0,1,0,1,1,1,1,1,1,0,1,1,1,1,1,0,1,0,0,0,0,0,1,0,1,1,0,0,0,0,0,1,0,0,1,0,0,1,1,1,1,0,0,1,0,0,0,0,1,0,1,1,1,1,1,0,0,0,0,1,0,1,0,0,1,1,1,1,0,1,0,1,1;1,1,1,1,1,1,1,1,0,0,1,1,1,1,0,0,0,0,1,0,0,1,0,0,0,0,1,1,0,0,1,0,0,1,1,1,1,1,0,0,0,0,1,0,1,1,0,0,0,0,0,1,0,0,1,0,0,0,1,1,1,1,1,0,0,0,0,1,0,1,0,1,1,1,0,0,1,1,1,0,1,1,0,0,1,1,0,1,0,0,0,0,0,0,0,1,0,0,1,0;1,0,1,0,0,1,1,0,0,1,1,1,0,1,1,0,0,1,1,1,1,1,0,0,1,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,1,1,1,0,1,0,0,1,1,0,1,1,1,1,0,0,0,1,0,1,1,0,1,0,1,1,1,0,0,1,1,0,1,0,1,0,0,1,1,1,1,1,0,0,0,1,1,1,0,1,1,1,1,0,1,0,1,1,1,1;1,0,0,1,0,1,0,0,1,0,1,0,1,1,0,1,0,0,0,1,1,1,0,1,1,1,1,1,0,1,1,1,1,1,1,1,0,0,0,1,0,0,1,0,1,1,0,0,0,1,1,0,1,1,0,0,1,0,1,1,0,1,1,0,0,1,1,0,0,0,0,0,0,1,1,1,0,0,0,1,1,1,1,0,1,1,1,0,0,1,0,1,1,1,1,0,0,1,0,0;1,1,0,1,1,0,1,1,0,1,1,1,1,1,1,0,0,1,0,0,0,0,0,0,0,1,0,1,0,1,0,1,1,0,0,0,0,1,1,0,0,1,1,1,0,1,1,1,1,0,0,0,0,1,1,0,0,1,0,0,0,0,0,1,1,0,1,1,1,1,0,1,1,0,1,0,1,0,0,0,0,0,1,0,0,1,0,1,0,1,1,1,0,1,0,1,1,1,0,0;0,0,0,0,0,1,0,1,1,1,0,1,0,1,1,0,0,0,1,1,1,0,1,1,1,1,0,0,1,1,1,1,1,0,1,0,1,0,0,0,1,1,0,1,0,1,1,1,1,1,0,0,1,1,0,1,0,1,0,0,1,0,1,1,1,1,1,1,0,0,0,0,1,0,0,0,1,0,1,1,0,1,1,0,0,1,0,0,1,0,1,0,1,1,0,0,0,0,0,0;1,0,1,0,0,1,1,1,0,1,0,0,1,0,0,1,0,0,1,1,0,0,0,1,0,1,1,0,1,0,1,0,1,0,0,0,0,1,0,1,1,1,0,1,0,0,0,1,0,0,1,1,0,0,1,0,1,0,1,0,1,0,1,1,1,1,0,0,1,1,0,0,0,1,0,0,0,1,0,1,1,1,0,0,1,1,1,1,0,0,0,0,1,0,1,1,1,0,0,1;1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,1,0,0,0,1,1,1,0,1,1,0,1,1,0,1,1,0,0,0,0,1,1,1,1,1,1,0,0,1,1,0,0,0,0,1,0,1,1,0,1,1,0,1,1,0,1,1,1,0,1,1,0,1,1,0,1,1,1,1,1,0,1,0,1,0,0,1,1,0,1,1,0,0,1,0,1,0,1,0;1,1,0,0,1,1,1,1,1,1,0,1,0,1,0,0,0,1,0,0,1,0,1,0,0,1,1,0,1,0,0,1,1,0,0,1,0,0,1,1,0,1,0,0,0,1,0,1,0,0,0,1,0,0,0,0,1,0,1,0,1,0,1,0,0,0,0,1,0,1,1,1,0,0,1,1,1,0,0,0,1,1,1,0,1,1,0,0,0,1,0,1,1,0,1,1,0,0,0,1;1,0,1,0,1,0,0,0,0,0,1,1,0,0,1,1,1,0,1,0,0,0,1,1,1,1,0,0,0,0,0,1,0,1,1,0,0,0,0,0,1,0,0,1,1,0,0,1,1,0,0,0,1,1,1,0,0,1,1,0,1,1,1,0,1,1,0,0,0,1,0,1,1,0,1,0,1,0,1,1,1,1,0,0,0,1,0,1,1,1,1,1,0,0,0,1,1,0,0,0;1,0,1,0,0,1,0,0,0,1,0,1,1,0,1,1,1,0,0,1,0,1,0,1,1,0,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,0,1,1,0,0,0,1,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,1,1,1,1,1,0,1,1,0,0,1,0,1,1,1,0,0,0,1,1,0,1,1,0,1,1,1,0,1,0,1,0,1,0,1,1;0,0,1,1,1,0,0,1,0,1,0,1,1,0,1,0,0,1,0,0,1,1,0,0,1,0,1,1,0,0,1,0,1,1,1,0,1,0,0,0,1,1,1,1,0,0,0,0,1,0,1,1,0,0,0,0,1,0,0,0,0,0,1,1,0,0,1,0,0,0,0,1,0,0,1,1,0,1,0,0,1,0,0,0,1,1,1,1,1,1,1,0,0,1,0,0,1,0,1,0;0,0,0,1,0,1,1,1,1,0,1,1,1,0,0,0,1,0,0,1,1,1,1,0,0,1,0,1,1,1,1,1,0,0,0,0,0,0,1,0,1,0,1,1,0,1,0,0,1,0,1,0,0,0,0,0,0,1,1,1,1,0,1,0,1,0,0,0,1,0,1,0,1,1,0,1,1,0,0,0,0,0,1,0,1,0,0,0,0,1,1,1,1,0,0,0,1,1,1,1;0,0,1,0,0,1,0,1,1,1,0,0,0,0,1,0,0,1,1,0,0,1,1,1,0,1,0,0,0,0,1,1,1,0,1,0,1,0,0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,1,0,0,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,0,1,0,1,1,0,1,1,0,1,0,1,1,1,1,1,0,0,1,0,0,1,1,0,1,1,0,1,1;0,0,0,0,1,0,0,0,1,1,0,0,1,0,1,1,0,0,1,1,0,0,1,1,1,0,1,0,0,0,1,0,1,1,1,0,1,0,0,1,0,1,1,1,1,0,1,0,0,1,0,0,1,0,1,0,1,0,0,0,1,0,1,0,0,0,0,1,1,0,1,1,0,0,0,1,1,0,0,1,1,1,0,0,0,0,1,0,1,1,0,1,1,0,0,0,1,1,0,0;1,1,1,0,1,0,1,1,0,0,0,1,1,0,1,0,0,0,1,1,1,0,0,1,1,1,0,1,1,1,1,0,1,1,0,0,0,0,0,0,0,0,0,1,0,1,1,1,0,1,0,0,1,0,0,1,0,0,1,1,0,0,0,1,1,0,0,0,1,0,0,0,1,0,0,0,0,1,1,0,0,0,1,1,0,0,1,0,0,0,1,0,0,1,0,0,0,1,1,0;0,1,1,0,1,0,0,1,0,0,0,0,1,1,1,1,1,1,1,0,0,1,1,0,1,1,1,0,0,0,0,1,1,1,0,1,0,0,0,1,1,0,1,0,0,1,0,1,1,0,0,0,0,1,1,0,1,1,0,1,1,0,0,0,1,0,1,1,1,1,1,0,0,1,0,1,0,1,0,1,0,0,0,0,0,0,0,0,1,0,0,1,1,1,0,1,0,1,1,1;1,0,1,0,0,1,0,0,0,1,1,0,1,0,0,1,1,1,0,1,1,0,0,1,0,1,1,1,0,1,0,0,1,0,1,1,1,1,0,1,0,1,0,0,1,0,1,1,1,1,1,1,0,0,1,1,1,0,0,0,1,0,1,0,0,0,1,1,0,1,1,1,1,0,0,0,0,1,0,1,0,1,1,0,1,1,1,0,1,0,0,0,0,0,1,1,1,1,1,1;0,1,1,0,1,0,1,1,0,0,1,0,1,1,0,0,1,0,0,1,1,1,0,0,1,0,1,1,1,0,0,1,0,1,1,1,1,0,1,1,1,0,0,0,1,0,0,0,0,1,0,1,0,1,0,1,1,1,1,1,1,1,1,1,0,1,0,1,0,0,0,0,1,1,0,1,1,0,1,0,0,0,0,0,0,0,1,1,0,1,0,1,1,1,1,1,1,0,1,1;0,1,0,0,0,0,1,0,0,1,0,0,0,0,1,1,0,1,0,1,0,1,0,0,1,0,0,1,1,1,0,1,1,0,1,0,1,1,0,1,1,1,0,1,1,0,0,1,0,1,0,1,0,0,1,0,0,0,0,0,1,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,1,0,0,0,0,0,1,0,1,1,1,0,1,0,1,1,1,0;1,0,1,0,0,1,0,0,1,0,0,0,1,1,1,0,1,0,0,1,0,1,0,0,1,0,1,0,1,1,0,0,1,1,1,1,1,1,0,1,0,1,1,1,0,0,0,0,0,0,1,1,0,1,1,1,0,0,1,0,1,1,0,1,0,1,1,0,0,0,0,1,1,1,1,0,0,1,0,1,0,1,0,0,0,0,0,0,0,1,0,1,1,0,0,0,1,1,1,1;1,1,0,1,0,0,1,0,0,1,1,0,1,0,1,1,1,0,0,1,1,1,1,1,1,0,0,0,0,0,1,0,1,1,0,1,0,1,0,0,1,0,0,1,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,1,1,0,0,0,0,1,0,0,0,1,1,0,0,0,1,1,1,0,1,0,1,1,1,0,1,1,0,0,1,1,0,0,1,1,0,0,1,0;1,0,1,1,1,1,1,0,0,1,0,1,1,1,1,0,0,1,1,1,0,1,1,0,0,1,0,1,1,0,0,1,1,1,1,0,0,0,0,1,1,1,0,1,1,0,1,0,1,0,0,0,0,0,1,0,0,1,1,1,0,1,0,1,0,1,1,1,1,0,1,1,0,0,1,0,0,1,0,0,1,1,0,1,0,0,0,1,0,1,0,0,1,0,1,0,0,1,0,1;1,1,0,1,0,1,0,0,1,1,0,0,1,1,1,1,1,1,0,1,1,0,1,1,1,1,1,0,1,1,1,1,1,1,1,1,0,1,1,0,1,0,1,0,1,0,1,1,0,0,1,1,1,0,0,1,0,1,0,0,1,0,0,1,1,1,1,1,0,1,1,1,1,0,0,0,1,0,1,0,0,0,0,1,1,1,0,1,0,1,1,1,0,0,1,0,0,1,1,1;0,1,0,1,0,0,1,0,0,1,1,1,1,0,0,0,1,0,1,1,1,0,0,1,1,1,0,1,0,1,1,1,1,1,0,1,0,1,1,1,0,1,0,1,0,1,0,0,1,1,0,1,1,0,0,1,1,0,1,0,0,1,0,0,0,1,1,1,0,0,0,1,1,0,1,1,0,0,0,1,1,1,0,0,0,0,0,1,0,0,0,1,1,0,1,0,1,0,1,0;0,0,0,0,0,0,1,0,1,1,1,0,1,0,1,0,0,0,1,1,1,0,1,1,0,0,1,1,1,1,0,1,1,0,1,1,1,0,1,1,0,1,0,0,0,1,1,1,0,1,0,0,1,0,0,0,0,1,0,0,0,1,1,1,1,0,0,0,0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,1,1,1,0,1,0,1,0,0,0,1,0,0,0,1,0,1;1,0,0,1,1,0,1,1,1,1,1,1,1,0,0,0,1,1,0,1,0,0,0,0,0,1,1,1,0,1,1,0,1,1,1,1,0,1,0,0,1,1,1,1,1,0,0,1,1,1,0,1,1,0,0,1,0,1,1,1,0,1,1,0,0,0,1,1,1,1,0,1,1,1,0,0,0,1,0,0,0,1,0,1,0,1,0,1,0,1,1,1,1,0,1,0,1,0,0,0;0,0,1,1,1,0,0,1,1,0,1,0,0,0,1,0,1,0,0,0,1,0,1,1,0,0,1,1,1,1,0,0,0,0,1,0,1,1,0,0,0,1,1,0,0,0,0,0,1,1,1,1,1,0,1,1,0,1,0,1,0,0,1,0,0,1,0,1,0,0,0,1,0,0,1,0,0,1,1,1,0,0,1,0,1,0,1,0,1,1,1,1,1,0,1,1,0,1,0,1;1,1,0,0,0,0,0,0,1,1,1,0,0,1,0,1,0,0,0,1,0,0,0,0,0,0,1,0,1,1,1,0,1,1,0,1,1,1,0,1,1,1,1,0,1,0,1,1,1,0,1,1,1,0,1,0,1,0,0,0,0,1,1,0,1,0,1,1,0,0,1,1,1,1,0,1,0,1,1,1,1,0,1,1,1,0,1,1,1,1,0,1,0,0,0,0,0,0,0,0;1,1,1,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,1,0,1,0,0,0,0,0,1,0,0,0,0,1,1,1,0,0,0,1,1,0,1,1,1,1,0,0,1,1,0,0,0,1,0,0,1,1,0,0,0,1,1,0,1,0,1,1,0,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,1,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1;1,0,1,0,0,1,0,1,1,1,0,0,1,0,0,1,0,1,0,1,0,0,0,1,0,1,1,1,1,1,0,1,0,1,1,0,0,1,1,1,0,1,1,1,0,0,0,1,1,0,1,1,1,0,1,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,1,0,0,0,1,0,1,0,1,1,1,1,0,0,0,0,0,0,0,0,1,0,1,1,1,1,1,1,0,1;1,1,1,0,1,1,0,1,0,1,0,1,0,0,1,1,0,0,1,1,1,1,0,0,0,1,0,1,1,0,1,1,1,0,0,1,0,1,0,0,1,1,0,0,0,1,1,0,0,1,1,0,0,0,0,1,1,0,0,1,1,1,1,1,1,1,1,1,0,1,0,0,1,0,1,0,1,0,1,1,1,0,0,0,1,0,1,0,0,1,0,0,1,0,0,1,0,1,0,1;1,0,1,0,0,0,0,0,1,1,0,1,0,1,1,1,1,1,0,1,1,0,0,1,0,0,1,0,1,1,0,1,0,1,1,1,1,1,1,1,1,1,0,0,0,1,1,1,0,1,1,0,0,1,0,0,0,0,0,0,1,1,0,0,1,1,0,0,0,1,1,0,1,1,1,1,1,0,0,0,0,0,1,0,1,1,0,0,0,1,0,0,1,1,0,1,1,0,1,1;0,1,1,1,1,0,0,0,0,0,1,1,1,1,0,0,1,0,0,0,1,1,0,1,0,1,0,0,0,1,0,0,1,0,0,1,1,1,1,1,0,0,1,0,0,1,0,0,1,1,1,0,0,1,1,1,1,0,1,1,1,0,1,1,1,0,0,1,0,0,0,1,1,0,1,1,1,0,0,1,1,1,1,0,0,1,1,0,0,0,1,1,1,1,0,1,1,1,1,1;1,1,1,1,1,0,1,0,0,1,0,0,0,1,1,1,1,0,1,1,1,1,1,1,1,0,0,0,1,1,1,1,0,1,0,1,0,0,1,1,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,1,1,1,1,0,1,0,1,0,0,0,1,1,1,0,0,0,1,1,1,0,0,1,0,1,0,1,0,0,0,1,0,1,1,0,0,0,0,0,1,1,1,1,0,1;0,0,0,0,1,1,1,1,1,0,1,1,1,0,0,0,1,0,1,1,0,0,0,1,0,0,1,1,1,0,1,1,1,0,0,1,0,1,1,0,0,0,0,0,1,1,1,0,1,1,0,1,1,1,1,1,1,1,0,1,1,1,1,1,0,0,1,0,0,1,1,1,1,1,0,1,0,1,1,1,0,1,0,1,0,0,0,0,0,0,1,0,0,1,1,0,0,0,1,1;0,0,0,0,1,0,1,1,0,0,1,0,1,1,1,0,1,1,0,0,0,1,0,0,1,1,0,0,0,0,1,0,0,1,1,0,0,0,0,0,1,1,1,0,1,1,1,0,1,1,0,1,1,1,1,1,0,1,0,0,0,1,1,0,0,0,1,1,1,0,0,1,0,1,0,1,1,0,1,1,0,0,0,0,0,0,0,1,1,1,0,0,1,0,1,1,1,1,1,1;0,0,1,1,0,0,0,1,0,0,0,0,0,1,1,0,1,0,0,0,0,0,1,1,1,0,1,0,0,0,0,1,1,0,1,0,0,1,0,0,1,1,0,0,1,1,1,0,1,1,1,1,0,0,0,1,1,1,1,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,1,0,0,1,1,0,1,1,1,1,1,1,0,1,0,1,1,1,1;1,1,1,1,0,0,1,1,0,0,0,1,0,1,1,1,0,1,1,0,0,0,1,0,1,1,1,0,1,0,0,0,1,0,1,1,0,1,1,1,0,1,0,0,0,0,0,1,0,0,0,0,0,0,1,0,1,0,0,0,1,0,0,0,1,1,1,0,0,0,0,1,1,0,0,0,1,1,0,0,0,1,0,1,1,0,1,1,1,0,1,0,0,1,0,0,0,0,1,0;1,0,0,0,1,1,1,0,0,0,0,1,0,1,1,0,0,0,1,1,1,1,0,0,0,1,1,0,0,0,0,1,0,1,0,1,1,1,1,1,0,0,1,1,1,1,1,0,1,0,1,0,1,1,0,0,0,0,1,0,1,0,1,1,1,1,0,0,0,1,1,1,1,0,1,1,1,1,0,1,1,1,1,1,1,0,0,0,0,0,0,0,1,0,0,0,1,1,1,0;0,1,1,0,0,1,0,0,1,1,0,0,1,0,1,0,1,0,0,0,0,0,0,1,1,0,1,1,1,0,0,0,0,1,1,1,1,0,0,0,1,1,1,0,1,1,1,0,0,1,0,0,1,1,0,1,0,0,0,1,0,0,0,0,0,1,0,1,1,0,1,0,0,1,1,0,0,0,0,1,0,0,0,1,1,0,0,0,1,1,0,0,0,0,1,1,1,0,0,1;1,0,1,0,1,1,0,0,1,0,0,1,1,0,0,1,1,0,0,1,1,1,0,0,0,0,1,0,0,1,0,0,1,0,0,0,1,1,0,1,1,1,1,1,0,0,1,0,1,0,1,0,0,1,1,1,1,0,1,1,0,0,0,1,1,1,0,0,0,0,1,1,1,0,1,0,1,1,1,0,0,1,1,0,0,1,1,0,0,1,0,0,1,1,1,1,1,1,1,1;0,0,1,0,0,0,0,0,0,1,1,1,0,0,0,1,0,1,0,1,1,0,0,0,0,0,1,1,1,1,0,0,1,1,0,1,1,1,0,1,0,0,0,0,1,1,1,0,0,0,0,1,0,0,0,0,0,1,1,1,0,0,0,1,1,0,0,0,0,0,0,0,0,1,1,1,0,1,0,1,0,0,1,1,1,0,0,0,1,1,1,0,0,0,1,0,1,0,0,1;1,1,0,1,0,1,0,1,1,1,0,1,1,0,1,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,1,0,0,1,1,0,0,1,0,1,0,0,0,1,1,0,1,0,0,0,1,1,0,1,1,1,0,0,0,0,0,0,0,1,1,1,1,0,1,0,1,1,0,1,1,0,0,0,1,0,1,0,0,1,1;1,0,0,0,1,0,1,0,0,0,0,1,1,1,1,0,0,0,1,0,0,0,1,0,0,1,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,0,1,1,1,0,0,1,1,0,0,1,0,1,0,1,0,1,0,1,1,1,0,0,1,0,1,0,1,1,1,1,0,1,0,1,0,1,1,1,0,0,1,1,0,1,1,1,0,0,1,1,1,0;0,0,0,0,0,1,0,0,1,0,1,0,0,1,0,1,0,0,1,0,0,0,0,1,0,1,1,0,1,1,0,1,0,0,0,0,1,1,0,1,0,0,1,0,1,1,0,1,0,0,1,0,1,1,1,1,0,0,0,0,0,1,0,0,1,0,1,0,0,0,1,1,1,0,0,1,1,0,1,0,1,0,0,0,1,0,0,0,1,0,1,1,1,1,0,0,0,0,0,0;1,1,1,0,0,1,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,1,1,0,1,0,0,1,1,0,1,1,0,1,0,1,0,1,1,1,1,1,0,0,1,1,0,0,0,1,1,0,1,1,0,1,0,0,0,1,0,0,0,1,0,0,1,1,1,1,1,1,1,0,1,0,0,0,0,0,0,1,1,0,0,1,1,0,1,0,1,1,0,0,0;0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,1,0,1,0,0,1,0,1,1,0,1,1,1,0,0,0,0,0,1,0,0,0,1,1,0,1,0,1,1,1,0,1,1,0,0,1,0,0,0,0,0,1,1,1,0,1,1,1,1,0,0,1,0,0,0,0,1,1,0,1,0,1,1,1,0,1,0,0,1,0,1,1,1,1,1,0,1,0,1,1,1,0,0,1,0;1,1,0,1,1,0,0,1,1,0,0,1,0,1,1,0,1,0,1,1,0,1,1,0,0,1,0,1,0,0,0,1,1,0,1,1,1,0,0,0,0,0,0,1,1,1,1,0,0,0,0,1,0,1,0,1,1,1,1,1,0,0,1,0,0,0,1,0,1,0,0,0,0,0,1,1,1,1,1,1,1,1,0,1,0,1,1,1,0,0,0,1,1,1,0,0,0,0,0,1;0,0,1,0,0,0,1,1,1,1,1,0,1,0,0,1,1,1,1,0,0,1,1,0,1,0,0,1,0,1,0,1,0,1,0,1,0,0,0,0,0,0,1,1,0,0,1,0,1,0,1,1,1,0,0,1,1,1,1,0,1,0,1,0,0,1,0,0,0,1,0,1,1,0,1,0,0,1,1,0,0,1,1,0,1,1,1,1,1,0,1,0,1,0,1,0,0,0,1,1;0,0,1,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1,1,0,1,1,0,1,0,0,0,1,0,0,0,1,1,0,0,0,1,0,1,0,1,0,0,0,0,1,1,1,1,1,0,0,0,1,0,1,1,1,1,0,0,1,1,1,1,0,0,0,1,1,1,1,1,0,1,1,1,0,1,0,0,1,1,1,0,0,0,1,0,1,1,1,1,1,1,1;1,1,0,0,0,0,1,0,0,0,1,1,0,0,1,1,1,1,1,0,0,1,1,1,0,1,1,1,1,1,1,0,1,0,0,0,0,0,1,1,1,1,1,1,1,0,1,1,1,0,0,0,0,0,0,1,1,0,1,1,1,1,1,0,1,0,1,0,0,1,1,0,1,1,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,1,0,0,1,1,1,1,1,0;1,0,0,1,0,1,0,1,0,1,1,0,1,0,0,0,1,0,1,0,0,0,0,0,0,0,0,1,0,1,1,1,0,1,1,1,0,1,1,0,1,1,0,0,1,1,0,0,0,0,0,0,1,1,1,0,1,0,0,1,1,1,0,0,1,1,1,1,1,0,0,1,0,0,0,1,0,1,1,1,0,0,1,0,0,1,1,1,0,1,1,1,0,0,0,0,0,1,0,1;1,0,1,1,1,1,0,0,0,0,1,1,1,0,1,1,0,1,1,0,1,1,1,1,0,0,1,1,1,0,1,0,0,0,1,1,1,1,0,1,1,0,1,1,1,1,0,0,1,0,0,0,0,0,0,0,1,1,1,1,1,0,1,1,1,1,1,0,0,0,0,1,1,1,1,1,0,0,1,0,1,0,0,0,1,0,0,1,0,0,0,0,1,1,1,0,0,1,1,0;1,0,1,0,0,0,0,0,1,0,0,0,0,1,1,1,1,0,0,0,1,0,1,0,1,0,0,1,1,1,0,1,1,0,1,0,0,0,1,0,1,0,1,0,1,0,1,0,1,0,1,1,0,1,0,0,1,0,0,0,0,0,1,1,0,1,1,1,0,1,0,1,0,0,0,1,0,1,1,1,1,0,1,0,1,0,0,1,1,0,1,0,0,1,1,1,0,1,1,0;1,0,0,1,0,0,0,1,1,1,0,1,0,1,1,1,1,0,1,1,0,1,0,0,1,1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,1,1,1,0,0,0,1,1,1,0,1,1,0,0,1,1,0,0,0,0,1,1,1,0,1,1,1,1,0,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,1,1,0,0,1,1,1,1,0,1,0,0,0,1,0,1;0,0,0,1,1,1,0,1,1,1,0,1,1,0,1,1,0,0,1,1,0,0,0,0,0,0,0,1,0,1,0,1,1,1,0,0,1,0,1,1,1,1,0,0,0,0,0,1,1,1,1,0,1,1,0,0,0,0,1,1,0,1,1,1,1,1,0,1,0,0,1,0,1,0,1,0,0,1,1,1,0,0,1,1,0,1,0,0,1,0,0,1,1,0,0,0,1,0,0,1;1,0,0,1,0,1,0,0,0,1,0,1,1,1,1,0,1,0,0,1,1,0,0,0,0,1,1,0,1,1,0,1,1,1,0,1,0,1,1,0,1,0,0,1,1,1,0,1,0,0,0,0,1,1,1,0,1,1,0,1,1,1,1,1,1,0,1,1,1,0,1,0,1,1,1,1,1,0,0,1,0,0,1,1,1,0,1,0,1,1,0,1,1,1,0,0,1,0,0,0;0,1,0,0,1,0,0,0,1,1,1,0,0,1,1,0,1,1,0,1,0,0,0,1,0,1,1,1,0,0,1,1,1,1,0,1,1,1,0,0,1,0,1,1,0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,1,0,1,0,1,1,1,1,1,1,1,0,0,1,1,1,0,1,1,1,1,1,0,1,0,1,1,0,1,0,0,1,0,1,1,1,1,0,0,0,1;0,1,0,1,1,1,0,0,0,0,0,0,0,1,0,1,1,0,0,1,0,1,0,1,1,1,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,0,1,0,0,1,0,1,0,1,0,1,0,0,0,0,1,1,1,1,0,0,0,1,0,0,0,1,0,0,1,1,0,0,0,0,1,0,0,1,1,0,1,1,0,1,1,1,0,1;0,0,1,0,1,1,0,1,0,0,1,1,0,1,0,1,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,1,0,0,1,0,1,1,0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,1,0,1,0,0,1,0,0,0,1,1,1,1,0,0,1,1,0,0,0,1,1,0,1,1,1,0,1,0,1,1,1,0,1,0,0,1,1,0,1,1,0;0,1,1,0,0,0,1,0,1,1,0,1,0,0,0,0,1,1,0,1,0,1,0,1,0,1,1,0,0,0,0,1,1,0,0,0,0,1,0,1,0,1,0,0,1,0,0,0,1,1,1,0,1,1,1,0,0,0,0,0,1,0,0,0,1,1,1,0,0,1,1,1,1,0,0,0,1,0,0,0,0,0,1,1,0,0,0,1,0,1,0,0,0,1,1,0,1,1,0,1;0,0,0,1,0,0,1,0,1,0,1,0,0,1,0,0,1,1,1,1,1,0,0,1,0,0,1,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0,1,0,1,1,0,1,1,0,1,0,0,0,1,1,1,0,1,0,0,1,1,1,0,0,0,0,0,0,1,1,0,0,0,0,1,1,0,1,0,1,1,1,0,1,0,0,0,0,0,1,0,0,0,1,1,1,0,1;1,0,0,1,1,0,1,0,1,0,1,1,0,1,1,0,0,0,1,0,0,1,1,0,1,0,1,1,0,1,1,0,1,1,0,1,0,1,0,0,1,1,1,1,1,0,0,1,1,0,1,0,0,1,1,1,1,0,1,1,1,0,1,0,0,1,1,1,0,0,1,0,1,0,0,0,1,1,1,0,1,0,1,1,0,0,0,0,1,1,1,0,0,1,1,0,0,0,0,0;1,1,0,0,1,0,1,1,1,0,1,0,1,0,0,1,1,0,0,0,0,1,0,0,0,1,0,1,0,1,0,0,0,0,1,1,0,1,0,0,0,1,0,1,1,1,0,0,0,1,0,1,0,0,0,1,0,0,0,1,1,0,1,0,0,0,1,1,1,1,0,0,0,1,1,0,0,1,1,1,1,1,0,1,0,0,0,1,0,1,0,1,0,0,0,1,1,1,0,0;1,1,0,0,0,0,1,0,1,0,0,1,1,1,0,0,1,1,1,1,1,0,1,0,0,0,0,0,0,1,0,1,0,1,0,0,1,0,1,1,1,1,1,1,0,0,0,0,1,1,1,1,0,1,0,1,1,1,1,1,1,0,1,0,0,1,1,1,0,1,0,0,0,1,1,0,1,1,1,0,0,1,0,1,1,1,0,0,1,0,0,1,1,1,1,0,1,0,1,0;0,0,0,1,1,1,1,1,1,1,0,0,1,0,0,0,0,1,0,0,1,1,1,1,0,1,0,1,0,0,0,0,0,1,0,0,0,1,0,0,0,1,1,0,1,1,0,0,1,0,0,1,0,1,1,1,0,1,0,1,0,1,1,1,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,1,1,0,0,1,1,0,1,0,0,1,0,1,0,0,1,0,1,1,0,0;0,0,1,0,1,0,0,1,1,0,1,0,0,1,1,0,1,1,1,1,0,1,0,1,0,0,0,1,1,0,1,0,1,0,0,0,0,0,0,1,1,1,1,0,0,1,1,1,1,0,1,0,0,1,1,1,1,1,0,1,0,0,0,0,0,0,1,1,0,0,1,1,1,0,1,0,1,0,1,0,0,0,1,0,1,1,1,1,1,1,0,0,0,1,0,0,1,0,1,0;1,1,1,0,1,0,1,1,1,1,1,1,0,0,0,1,1,0,0,1,1,0,1,0,1,1,1,0,1,1,1,1,0,0,0,1,1,1,0,0,0,0,0,1,1,0,0,1,1,0,1,1,0,1,0,1,1,1,1,0,0,1,0,1,1,1,0,1,1,0,0,1,1,1,1,0,0,1,0,0,1,1,0,0,1,1,0,1,1,1,1,0,1,1,1,0,0,0,1,1;0,1,1,0,1,0,1,1,0,1,1,1,0,0,1,0,1,0,1,1,0,0,1,0,1,0,0,1,0,0,1,0,1,0,1,0,1,1,0,1,1,0,0,0,1,1,1,0,0,0,1,0,1,0,1,0,1,1,1,1,0,1,1,1,0,1,0,1,0,1,0,0,1,1,1,0,1,0,1,0,1,0,0,1,1,0,1,0,0,0,1,1,0,1,0,0,1,1,0,0;1,0,1,1,0,1,0,0,0,1,0,1,1,0,1,1,1,0,1,0,0,0,0,1,0,1,1,0,1,1,0,0,0,1,0,0,1,1,0,1,1,0,1,1,1,1,0,0,1,1,0,1,1,1,0,1,0,1,0,1,0,1,0,1,1,1,1,1,0,1,0,1,0,1,0,1,0,0,0,1,1,1,1,0,0,1,1,1,1,1,0,1,0,1,1,1,0,0,0,0;1,0,0,0,0,1,1,0,0,1,1,1,1,0,0,1,1,1,1,0,1,0,1,1,0,0,0,0,1,0,1,1,0,1,0,0,0,1,1,1,1,0,1,0,0,0,1,0,1,0,0,0,1,0,1,0,1,1,0,1,1,0,1,1,0,0,0,1,1,0,0,0,1,1,0,1,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,0,1,1,1,1,1,0;0,0,0,1,1,1,0,1,0,1,1,1,1,0,1,1,0,1,1,0,0,0,0,1,0,0,1,0,0,1,0,1,0,1,0,1,0,0,0,1,0,0,1,1,1,0,1,1,1,0,1,0,1,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,1,1,0,1,0,1,1,0,0,1,0,1,1,1,1,0,1,0,0,1,1,1,1,0,0,0,0,0,0,0,1,0;0,1,0,0,0,0,0,1,1,0,0,0,1,1,1,0,1,1,0,1,0,1,1,0,1,0,1,0,1,0,1,0,0,0,0,0,1,1,0,0,0,1,1,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,1,1,1,1,0,1,1,1,1,1,0,1,1,1,1,0,0,0,1,0,0,1,1,1,1,1,1,0,0,0,0,0,1,1,0,1,0,1,0,0,0,0;0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,0,1,1,1,0,1,1,0,1,1,0,0,0,0,0,1,0,0,1,1,1,0,1,1,1,0,0,1,1,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,1,1,1,0,0,1,0,1,0,1,1,1,0,0,0,1,1,1,0,1,1,1,1,0,0,1,0;1,1,0,0,0,1,0,1,0,0,1,0,1,0,0,1,1,1,0,0,1,1,1,0,0,0,1,0,0,0,1,0,1,0,1,0,1,1,1,0,1,1,0,0,0,0,1,1,1,1,0,1,0,1,1,0,0,0,1,0,0,0,1,1,1,0,1,1,0,0,0,0,0,0,1,1,1,1,1,0,1,1,1,1,1,1,1,1,0,1,1,0,0,1,0,1,1,1,0,0;0,1,1,1,0,0,0,0,1,0,1,1,1,1,1,1,0,1,1,1,1,0,1,0,0,0,1,0,0,0,0,0,1,0,1,1,0,0,0,0,0,1,1,1,0,0,1,0,0,0,1,0,1,1,0,0,1,1,1,1,0,1,0,0,1,1,0,1,0,1,0,1,0,0,1,0,1,1,0,1,0,0,0,0,1,1,0,1,0,0,1,1,1,0,1,1,0,0,1,1;0,1,1,1,0,1,1,0,1,1,0,1,1,0,0,1,0,0,0,1,1,0,1,1,1,0,1,1,0,0,1,0,0,0,0,0,1,1,1,0,1,0,1,0,0,0,0,1,0,0,1,0,1,0,0,1,1,1,1,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,1,0,1,1,0,0,0,0,1,0,1,1,1,0,0,0,0,0,0,0,0,1,0,1;1,0,1,1,1,1,0,0,0,0,1,1,0,1,0,1,1,0,1,0,1,0,0,0,0,0,0,1,0,0,1,1,1,1,1,1,0,1,1,0,0,0,0,1,0,1,1,1,0,0,0,0,0,0,0,1,1,1,1,1,0,1,1,1,0,0,0,1,0,1,1,0,0,1,0,0,1,1,0,1,0,1,0,0,1,1,1,1,0,1,1,0,0,0,1,0,0,0,0,0;0,1,0,0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,1,1,1,0,0,1,0,1,1,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,1,1,1,0,1,0,1,1,1,1,0,1,0,1,0,0,0,0,1,1,1,1,0,0,1,0,0,1,0,1,0,1,1,0,1,0,1,0,1,0,0,1,0,1,1,1,1,0,1,0,0,0,1,0,0;0,1,0,1,0,0,0,1,0,0,0,1,1,1,0,0,0,1,1,1,1,1,1,1,0,0,0,1,0,1,0,1,1,0,1,1,1,1,0,0,1,1,0,0,0,1,1,0,0,1,1,1,1,1,0,0,1,0,0,0,0,1,0,0,1,0,1,0,1,1,1,0,1,1,0,1,1,1,0,1,0,1,0,1,1,0,0,1,1,1,1,1,1,0,1,0,1,0,0,1;0,1,1,0,1,0,1,1,0,0,0,1,0,1,1,0,1,0,1,1,1,1,0,0,1,0,0,0,1,0,1,0,1,0,0,1,1,0,1,1,0,0,1,0,1,0,1,1,0,0,0,1,0,0,1,1,0,0,1,0,0,1,0,1,1,0,0,1,1,0,0,0,1,0,0,0,0,1,1,0,0,1,1,1,1,1,0,1,1,1,1,1,0,1,1,1,0,1,0,0;1,0,0,1,1,1,0,0,1,1,0,1,1,1,0,0,1,1,1,0,0,1,0,1,0,1,0,1,1,1,1,0,1,1,0,1,1,1,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1,1,1,1,1,0,1,1,1,0,0,1,1,1,0,0,1,0,1,0,1,1,1,0,0,1,1,0,0,1,0,0,1,0,0,1,1,1,1,0,1,0,1,1,0,0,0;1,0,0,0,0,0,0,1,1,1,0,1,1,0,1,1,0,1,0,1,0,1,1,1,0,1,0,1,1,1,0,1,0,1,0,1,1,0,0,1,1,1,1,0,0,1,1,0,1,0,1,0,0,1,1,0,0,1,1,0,0,0,1,0,0,1,1,1,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,1,0,1,0,0,0,1,0,0,1,0,0,1,1,1,0,0;0,1,1,1,0,0,0,1,1,1,0,0,1,1,1,0,0,0,0,0,1,0,1,0,1,1,0,1,0,0,0,0,0,0,1,0,0,0,0,1,0,1,1,0,1,0,0,1,0,0,1,0,1,1,1,1,1,1,0,1,0,0,1,1,1,0,1,1,1,0,1,0,1,0,1,0,1,1,1,1,0,0,1,1,1,0,0,0,1,0,1,1,0,1,0,0,1,1,1,1;1,1,1,0,1,1,1,1,1,1,0,1,1,0,0,1,1,1,0,1,0,0,0,0,0,0,1,1,1,0,1,1,1,1,0,1,1,0,1,1,0,0,0,1,1,1,1,0,0,1,1,1,0,0,0,0,1,0,1,1,1,0,1,1,0,0,0,1,0,1,1,0,1,0,1,1,0,1,0,1,1,0,0,1,0,1,0,1,0,1,1,0,0,0,1,1,1,1,1,0;0,1,0,1,1,1,1,0,1,0,1,0,0,1,0,1,0,1,1,0,0,0,1,0,0,1,1,1,0,0,1,0,0,0,0,0,1,0,1,1,1,1,1,1,0,1,0,0,0,1,1,0,1,0,0,1,1,0,0,1,1,0,0,1,0,0,0,1,1,1,0,1,0,1,0,0,0,0,0,1,1,0,1,1,1,1,0,0,0,0,1,1,1,0,1,1,0,0,0,1;0,1,0,1,1,0,0,0,0,1,0,1,0,1,0,1,1,0,1,1,1,1,1,1,0,0,1,1,1,1,0,0,0,1,0,1,0,0,0,1,0,1,1,1,0,1,1,0,1,1,1,1,0,1,0,1,0,0,0,1,1,0,0,0,0,1,1,0,1,0,1,1,0,1,1,1,1,0,1,0,1,0,0,0,1,0,0,0,0,1,0,1,1,1,1,0,1,1,1,1;1,1,1,0,1,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,1,0,1,1,1,1,0,1,1,0,1,1,0,1,0,1,0,0,1,1,0,1,1,0,1,1,0,1,0,1,0,0,1,0,0,0,0,0,1,1,1,1,1,1,0,0,0,1,1,1,1,0,1,0,1,0,0,1,0,1,0,0,0,1,0,1,0,1,0,1,0,1,1,1,0,1,1,0,0;0,0,0,0,1,0,0,0,0,1,1,1,0,0,0,0,1,0,0,1,1,1,1,0,1,1,1,1,1,1,1,0,1,1,0,0,0,0,1,0,0,1,1,0,1,1,1,1,1,0,1,0,1,1,0,0,1,0,1,1,1,0,1,1,0,0,0,0,0,1,0,0,0,0,1,0,1,1,0,0,1,1,0,1,0,1,0,0,0,0,0,0,0,1,1,0,1,0,1,1;0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,1,0,1,0,1,0,1,1,0,0,1,1,1,0,1,0,1,1,0,1,0,1,0,1,1,1,1,1,1,1,1,1,0,0,1,1,1,1,0,0,0,0,1,1,1,0,1,0,0,1,1,0,1,1,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,1,0,0,0,1,0,1,1,0,1,1];

%Generate a sequence of averaging matrices {A(k)}
mask_sequence = zeros(n , n , t_end);

for i = 1:t_end
   [mask_sequence( : , : , i ) , average_mask] = mask_samples_arbitrary_net(topology , n , p);
end

A_sequence = (1/(p+1))*mask_sequence;
average_A = (1/(p+1))*average_mask;

% Checking if the expected value is row stochastic;
row_stochastic = sum(average_A , 2);

end

%% EVOLUTION OF THE STATE
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

% Eigenvalue analysis of the E{A(k)}
%average_A is positive and has only one eigenvalue with |.|=1
[V , D , W] = eig(average_A);
lambda = eig(average_A);
v = V(: , end);
w = W(: , end);
 
% Limit of (E{A})^k is:
average_A_infinity = v*w';

%plotting eigenvalues of E{A}
if plotting_eig
    E_A = 'Eigenvalues of $\mathcal{E}[A(k)]$';
    figure(20);
    circle(0 , 0 , 1); 
    hold on; grid on; title(E_A , 'interpreter' , 'latex'); axis([-1.26 1.26 -1 1]);
    plot(eig(average_A) , zeros(n,1), 'r x', 'MarkerSize' , 3);
    hold off;
end

if lap
    % Computing the Laplacian matrix associated to the sequence of A(k)
    L_sequence = zeros(n , n , t_end);
    for i = 1:t_end
        L_sequence(: , : , i) = eye(n) - A_sequence(: , : , i);
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
end


%% IMPLEMENTING THE FEEDBACK CONTROL LAW
% number of coordinators
n_selfish = 5;

% Refererence chosen to be the mean of the standard agents
ref=mean(x_0(n_selfish+1:end))*ones(n_selfish , 1);

% Reference sequence (for plots)
ref_seq = mean(x_0(n_selfish+1:end)) * ones(t_end+1 , 1);

%Initializing Gains
Kp = 0.01*eye(n_selfish);
Ki = 0.001*eye(n_selfish);

%% STEP 1 P CONTROLLER
% Assumptions:
% - n_selfish coordinators IOTA NODES that are P-controlled
% - Expected Averaging Matrix
% - Global Selfish Agent,i.e. its control action is based on the whole network
% - Underlying network is a complete graph
% - The mean has to converge to a reference ref

if complete
    [x_k_average_P , y_k_average_P , average_A_P] = P_global(n , p , t_end , x_0 , n_selfish , ref , Kp);
end

%% STEP 1.1
%Same assumptions as in STEP 1, but now we introduce a saturation to the selfish node
%opinion so that it can't go higher than 1 or lower than 0

if complete
    [x_k_average_P_sat, y_k_average_P_sat, average_A_P] = P_global_saturation(n , p , t_end , x_0 , n_selfish , ref , Kp);
end

%% STEP 2 PI CONTROLLER
% Assumptions:
% - n_selfish coordinators IOTA NODES that are PI-CONTROLLED
% - Randomized Adjacency
% - Global Selfish Agent,i.e. its control action is based on the whole network
% - Underlying network is a complete graph
% - The mean has to converge to a reference ref

if complete
    [x_k_average_PI , y_k_average_PI , average_A_PI , average_A_PI_np] = PI_global(n , p , t_end , x_0 , n_selfish , ref , Kp , Ki);
end

% Checking if average_A_PI is Hurwitz, i.e. convergent
if plotting_eig
    eigen_average_A_PI_np = eig(average_A_PI_np);
    eigen_average_A_PI = eig(average_A_PI);
    if plotting_eig
        figure(21);
        circle(0 , 0 , 1); 
        hold on; grid on; title( 'Eigenvalues of closed loop system using Expectation'); axis([-1.26 1.26 -1 1]);
        plot(eig(average_A_PI) , 'r.', 'MarkerSize' , 1);
        hold off;
    end
end


%% Checking values of the gains that make the system in STEP 2 unstable
% if kp = 10 then ki can range between 0 and 10.22, before we go unstable
% if kp = 1 then ki can range between 0 and 0.8, before we go unstable
if plotting_eig
    figure(22);
    circle(0 , 0 , 1); 
    hold on; grid on; title( 'Checking values of k_I and k_P that make STEP 1 unstable'); axis([-1.26 1.26 -1 1]);
    for kp = 0:10
        for ki = 0:0.5:1
            average_A_cl_tuned = [average_A-B*kp*C , B*ki ; -C , 1];
            plot(eig(average_A_cl_tuned) , '.', 'MarkerSize' , 1);
        end
    end
    hold off;
end

%% STEP 2.1
%Same assumptions as in STEP 2, but now we introduce a saturation to the
%selfish nodes opinion so that it can't go higher than 1 or lower than 0

if complete
    [x_k_average_PI_sat, y_k_average_PI_sat, average_A_PI , average_A_PI_np] = PI_global_saturation(n , p , t_end , x_0 , n_selfish , ref , Kp , Ki);
end


%% STEP 3
% Assumptions:
% - n_selfish coordinators IOTA NODES that are P-CONTROLLED
% - Randomized Adjacency
% - Global Selfish Agents, i.e. its control action is based on the whole network
% - complete graph as the underlying network
% - The mean has to converge to a reference ref


[x_k_P , y_k_P , C] = P_rand(n , p , t_end , x_0 , n_selfish , ref , A_sequence , topology, complete, Kp);


%% STEP 3.1
% Same as 3 but with saturation on the opinions

if complete
    [x_k_P_sat , y_k_P_sat] = P_global_rand_sat(n , p , t_end , x_0 , n_selfish , ref , A_sequence , Kp);
end

%% STEP 4
% Assumptions:
% - n_selfish coordinators IOTA NODES that are PI-CONTROLLED
% - Randomized Adjacency used now instead of E{A}
% - Global Selfish Agents, i.e. its control action is based on the whole network
% - complete graph as the underlying network
% - The mean has to converge to a reference ref

if complete
    [x_k_PI , y_k_PI] = PI_global_rand(n , p , t_end , x_0 , n_selfish , ref , A_sequence , Kp , Ki);
end

%% STEP 4.1
% Same as 4 but with saturation on the opinions

if complete
    [x_k_PI_sat , y_k_PI_sat] = PI_global_rand_sat(n , p , t_end , x_0 , n_selfish , ref , A_sequence , Kp , Ki);
end

%% STEP 3
% Assumptions:
% - One selfish agent IOTA NODE
% - Randomized Adjacency
% - Myopic Selfish Agent C not strictly positive C = [c_1 , ... , c_n] only
%   m c_i's are non-zero, with m<n
% - complete graph as the underlying network
% - The error wrt the other nodes has to converge to a reference ref = 0

% Define the closed loop system state space representation (A , B , C , D)
B = [ones(n_selfish , 1) ; zeros(n - n_selfish , 1)];
C = [n-n_selfish , - ones(1 , n-n_selfish)];
Kp = 0.01;
Ki = 0.001;

% Defining new reference
ref = 0.001;

A_cl_sequence = zeros(n + n_selfish , n + n_selfish , t_end);
for k  = 1:t_end
    A_cl_sequence(: , : , k) = [A_sequence(: , : , k)-B*Kp*C , B*Ki ; -C , 1];
    x_k_PI(: , k+1) = A_cl_sequence(: , : , k) * x_k_PI(: , k) + [B*Kp ; 1]*ref;
end
figure(9) ; plot(0:1:t_end , x_k_PI(:,:) ,  'LineWidth' , 2); grid on; hold on; legend();
plot(0:1:t_end , x_k_PI(n+1 , :) ,  'LineWidth' , 2);title('Global, "Laplacian" reference'); hold off;
