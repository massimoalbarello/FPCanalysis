%% RANDOMIZED MODEL ON A COMPLETE GRAPH
% Defining the Adjacency matrix dimensions (n nodes)
n = 100;
% Defining the number of nodes to be sampled
p = 5;
%Fix simulation length
t_end = 100;       
%plotting variable
plotting = false;
plotting_eig = false;
lap = false;
complete = false;
not_complete = true;
first = true;

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
    average_G = digraph(average_A);
    
%     [pagerank,indices] = maxk(centrality(average_G,'pagerank'),n_selfish);
    

    %Check row stochasticity of A (satisfied)
    row_sums = sum(average_A,2);
    col_sums = sum(average_A,1);
    row_sum = sum(A_sequence,2);
    col_sum = sum(A_sequence,1);


%% GENERATING AN ARBITRARY TOPOLOGY (LOOK AT GRAPH_GEN FOR MORE DETAIL), THE NODES CAN SAMPLE FORM THE NEIGHBOURS ONLY, BUT ALWATS SAMPLE THEMSELVES
elseif not_complete
 
   if first
        % Variable to check if the requirements on the topology are met
        bad_topology = true;

        while bad_topology
            % Generating a suitable underlying graph
            [topology] = gen_graph2(n,p);

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


            if ((Summation==0) && (Irreducibility == 1))
                fprintf('Need more neighbours, Irreducible\n');
            elseif ((Summation==1) && (Irreducibility == 1))
                fprintf('Enough neighbours, Irreducible\n');
                topology = topology;
                bad_topology = false;
            elseif ((Summation==0) && (Irreducibility == 0))
                fprintf('Need more neighbours, Not Irreducible\n');
            elseif ((Summation==1) && (Irreducibility == 0))
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
        average_G = digraph(average_A);
   end
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


%% IMPLEMENTING THE FEEDBACK CONTROL LAW
% number of coordinators
n_selfish = 2;

% Refererence chosen to be the mean of the standard agents
ref=mean(x_0(n_selfish+1:end))*ones(n_selfish , 1);

% Reference sequence (for plots)
ref_seq = mean(x_0(n_selfish+1:end)) * ones(t_end+1 , 1);

%Initializing Gains
Kp =0.5*eye(n_selfish);
Ki = 0.01*eye(n_selfish);

%% DEFINING A TOPOLOGY THAT ALSO HAS A SINGLE MALICIOUS NODE (STUBBORN NODE). THE NETWORK WILL NOW HAVE N+1 NODES
n_stubborn = 1;
n_tot = n+n_stubborn;

if first
    topology_stubborn = zeros(n_tot , n_tot);
    new_connections = zeros( n-n_selfish ,n_stubborn);
    for i = 1:length(new_connections)
        new_connections( i , 1 ) = rand(1);
        if  new_connections(i , 1) > 0.9
            new_connections(i , 1) = 1;
        else
            new_connections(i , 1) = 0;
        end
    end

topology_stubborn = [topology(1:n_selfish , :), zeros(n_selfish , n_stubborn); topology(n_selfish+1:end , :) , new_connections ; zeros(n_stubborn , n), ones(n_stubborn)];
end
% Vector of initial conditions with stubborn agents
x_0_stubborn = zeros(n_tot , 1);
x_0_stubborn(1:n) = x_0;
% We want the malicious nodes to have an opinion that is the opposite of
% the majority of the network, thus they will all have opinion x = 0 and
% will never change it
x_0_stubborn(n+1:end , 1) = zeros(n_stubborn , 1);

if complete
    mask_sequence_stubborn = zeros(n + n_stubborn , n+n_stubborn , t_end);
    A_sequence_stubborn = zeros(n + n_stubborn , n+n_stubborn , t_end);
    for k = 1:t_end
        [mask_sequence_stubborn( : , : , k ) , average_mask_stubborn] = mask_samples_arbitrary_net_stubborn(topology , n , n_stubborn , p);
    end
    A_sequence_stubborn(1:n , : , :) = (1/(p+1))*mask_sequence_stubborn(1:n , : , :);
    A_sequence_stubborn(n+1:end , : , :) = (1/(n_stubborn))*mask_sequence_stubborn(n+1:end , : , :);
end
 

if not_complete && first
    mask_sequence_stubborn = zeros(n + n_stubborn , n+n_stubborn , t_end);
    A_sequence_stubborn = zeros(n + n_stubborn , n + n_stubborn , t_end);
    for k = 1:t_end
        [mask_sequence_stubborn( : , : , k ) , average_mask_stubborn] = mask_samples_arbitrary_net_stubborn(topology_stubborn , n , n_stubborn , p);
    end
    A_sequence_stubborn(1:n , : , :) = (1/(p+1))*mask_sequence_stubborn(1:n , : , :);
    A_sequence_stubborn(n+1:end , : , :) = (1/(n_stubborn))*mask_sequence_stubborn(n+1:end , : , :);
end


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
%Same assumptions as in STEP 1, but now we introduce a saturation to the
%opinion so that it can't go higher than 1 or lower than 0

if complete
    [x_k_average_P_sat, y_k_average_P_sat] = P_global_saturation(n , p , t_end , x_0 , n_selfish , ref , Kp);
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
% - Myopic Selfish Agents, i.e. limited visibility of the net
% - connected graph as the underlying network
% - The mean has to converge to a reference ref

[x_k_P , y_k_P ] = P_rand(n , p , t_end , x_0 , n_selfish , ref , A_sequence , topology, complete, Kp);

%% STEP 3.1
%Same assumptions as in STEP 3, but now we introduce a stubborn agent

[x_k_P_stubborn , y_k_P_stubborn] = P_rand_stubborn(n , p , t_end , x_0_stubborn , n_selfish , n_stubborn, ref , A_sequence_stubborn , topology_stubborn, complete, Kp);
    

%% STEP 3.2
% Same as 3 but with saturation on the opinions

%----SKIP----

%[x_k_P_sat , y_k_P_sat] = P_rand_sat(n , p , t_end , x_0 , n_selfish , ref , A_sequence , topology, complete, Kp);


%% STEP 3.3
% Same as STEP 3.2 but with stubborn agent disturbance

[x_k_P_sat_stubborn , y_k_P_sat_stubborn] = P_rand_sat_stubborn(n , p , t_end , x_0_stubborn , n_selfish , n_stubborn, ref , A_sequence_stubborn , topology_stubborn, complete, Kp);


%% STEP 4
% Assumptions:
% - n_selfish coordinators IOTA NODES that are PI-CONTROLLED
% - Randomized Adjacency used now instead of E{A}
% - Myopic Selfish Agents, i.e. its control action is based on the whole network
% - connected graph as the underlying network
% - The mean has to converge to a reference ref

[x_k_PI , y_k_PI] = PI_rand(n , p , t_end , x_0 , n_selfish , ref , A_sequence , topology , complete , Kp , Ki);

%% STEP 4.1
% Assumptions:
% - n_selfish = 1 coordinators IOTA NODES that are PI-CONTROLLED
% - Randomized Adjacency used now instead of E{A}
% - Myopic Selfish Agents, i.e. its control action is based on the whole network
% - connected graph as the underlying network
% - The mean has to converge to a reference ref
% - One stubborn agent

[x_k_PI_stubborn , y_k_PI_stubborn] = PI_rand_stubborn(n , p , t_end , x_0_stubborn , n_selfish , n_stubborn, ref , A_sequence_stubborn , topology_stubborn , complete , Kp , Ki);


%% STEP 4.2
% Same as STEP 4 but with saturation

[x_k_PI_sat , y_k_PI_sat] = PI_rand_sat(n , p , t_end , x_0 , n_selfish , ref , A_sequence , topology , complete, Kp , Ki);


%% STEP 4.3
% Same as 4.1 but with saturation
[x_k_PI_stub_sat , y_k_PI_stub_sat] = PI_rand_sat_stubborn(n , p , t_end , x_0_stubborn , n_selfish , n_stubborn, ref , A_sequence_stubborn , topology_stubborn , complete , Kp , Ki);

%% STEP 5
%Assigning coordinators according to centrality measures