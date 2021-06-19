function G = unweightedGraph(n)
%UNWEIGHTEDGRAPH Summary of this function goes here
    A = zeros(n);
    for i = 1:n
        for j = 1:n
            th = rand(1,1);
            %rand(n,n)
            if th > 0.5
                A(i,j) = 1;
                A(j,i) = 1;
            else
                A(i,j) = 0; 
                A(j,i) = 0;
            end
        end
    end
    G = graph(A);
end

