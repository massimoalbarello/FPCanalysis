function A_permuted = permute_matrix(centrality_vector,n_c,Adj_matrix)

n=size(Adj_matrix,1);
max_centrality = zeros(n_c,2); %max value in first column, index in the second

[max_centrality(:,1),max_centrality(:,2)] = maxk(centrality_vector,n_c);

P_matrix = eye(n);

for i=1:n_c
    row = zeros(n,1);
    row(max_centrality(i,2)) = 1;
    P_matrix(:,max_centrality(i,2)) = P_matrix(:,i);
    P_matrix(:,i)=row;
end    

A_permuted = P_matrix*Adj_matrix*P_matrix';

end