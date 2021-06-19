function G = simple_k_regularGraph(n,k)
%SIMPLE_K_REGULARGRAPH Summary of this function goes here
    % k = degree;
    matIter = 10;

    % check parameters
    if mod(n*k,2)==1   
        disp('input error: n*d must be even');
        A=[];
        return;
    end

    % graph adjacency matrix
    A=sparse(n,n);
    
    %a list of open half-edges
    U = repmat(1:n,1,k);
    
    edgesTested=0; 
    repetition=1;

    %continue until a proper graph is formed
    while ~isempty(U) && repetition < matIter

        edgesTested = edgesTested + 1;

        %print progress
        if mod(edgesTested, 5000)==0 
            fprintf('progress: edges=%k/%k\N', edgesTested, n*k);    
        end

        %chose at random 2 half edges
        i1 = ceil(rand*length(U));
        i2 = ceil(rand*length(U));
        v1 = U(i1);
        v2 = U(i2);

        %check that there are no loops nor parallel edges
        if (v1 == v2) || (A(v1,v2) == 1)

            %restart process if needed
            if (edgesTested == n*k)           
                repetition=repetition+1;            
                edgesTested = 0;
                U = repmat(1:n,1,k);
                A = sparse(n,n);
            end
        else
            %add edge to graph
            A(v1, v2)=1;
            A(v2, v1)=1;

            %remove used half-edges
            v = sort([i1,i2]);
            U = [U(1:v(1)-1), U(v(1)+1:v(2)-1), U(v(2)+1:end)];
        end
    end

    %check that A is indeed simple regular graph
    msg=isRegularGraph(A);
    if ~isempty(msg)    
        disp(msg);
    end
    
    G = graph(A);
end

function msg=isRegularGraph(G)
    %is G a simple d-regular graph the function returns []
    %otherwise it returns a message describing the problem in G

    msg=[];

    %check symmetry
    if (norm(G-G','fro')>0)
        msg=[msg,' is not symmetric, '];
    end

    %check parallel edged
    if (max(G(:))>1)
        msg=[msg,sprintf(' has %d parallel edges, ',length(find(G(:)>1)) )];
    end

    %check that d is d-regular
    d_vec=sum(G);
    if min(d_vec)<d_vec(1) || max(d_vec)>d_vec(1)
        msg=[msg,' not d-regular, '];
    end

    %check that g doesn't contain any loops
    if (norm(diag(G))>0)
        msg=[msg,sprintf(' has %d self loops, ',length(find(diag(G)>0)) )];
    end
end

