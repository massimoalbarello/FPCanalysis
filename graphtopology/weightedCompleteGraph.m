function G = weightedCompleteGraph(n)
%WEIGHTEDCOMPLETEGRAPH Summary of this function goes here
    A = zeros(n);
        for i = 1:n
            for j = 1:n
                th = rand(1,1);
                A(i,j) = th;
                A(j,i) = th;
            end
        end  
        G = graph(A);
end

