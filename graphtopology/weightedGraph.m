function G = weightedGraph(n)
%WEIGHTEDGRAPH Summary of this function goes here
 A = zeros(n);
    for i = 1:n
        for j = 1:n
            if i == j
                continue
            end
            th = rand(1,1);
            if th > 0.5
                A(i,j) = 1;
                A(j,i) = 1;
            else
                A(i,j) = 0; 
                A(j,i) = 0;
            end
        end
    end
    for i = 1:n
        for j = 1:n
            w = rand(1,1);
            A(i,j) = w*A(i,j);
            A(j,i) = w*A(j,i);
        end
    end
    G = graph(A);
end

