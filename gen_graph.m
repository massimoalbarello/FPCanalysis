function [topology] = gen_graph(n,p)
%Generate underlying topology like:
% - self loops
% - at least p+1 out neighbours, including himself
% - symmetric
% - connected
topology = zeros(n);
topology = randi(2,n,n) - 1;
topology = topology - tril(topology,-1) + triu(topology,1)'; % Building a binary simmetric matrix

% Introducing self loops
for i=1:n
    if topology(i,i) == 0
        topology(i,i) =1;
    end
end

%Checking irreducibility of topology binary matrix
Irreducibility = eye(n);
for i = 1:(n-1)
    Irreducibility = Irreducibility + topology^i;
end

% Checking the number of out neighbours including i itself, and irreducibility, if not good,
% we call the function again in recursive fashion
sum_row = sum(topology , 2);
if (sum_row < (p+1)*ones(n,1) & Irreducibility > zeros(n))
    fprintf('Need more neighbours, Irreducible\n');
    topology = gen_graph(n,p);
elseif (sum_row < (p+1)*ones(n,1) & Irreducibility <= zeros(n))
    fprintf('Need more neighbours, Not Irreducible\n');
    topology = gen_graph(n,p);
elseif (sum_row >= (p+1)*ones(n,1) & Irreducibility <= zeros(n))
    fprintf('Enough neighbours, Not Irreducible\n');
    topology = gen_graph(n,p);
elseif (sum_row >= (p+1)*ones(n,1) & Irreducibility > zeros(n))
    fprintf('Enough neighbours, Irreducible\n');
% else
%     fprintf('I forgot one case\n');
end

   
end

