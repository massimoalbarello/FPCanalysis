function graph_matrices = graphMatrices(G)
    %GRAPHMATRICES Summary of this function goes here

    % Description:  Given a graph, computes adjacency matrix, Laplacian, 
    %               incidence matrix, degree matrix
    % Output:       Data structure called graph_props that contains following properties of generated graph
    % 					graph_matrices.adjMat           : adjacency matrix
    % 					graph_matrices.degMat           : (out/in)degree matrix
    % 					graph_matrices.incMat           : incidence matrix
    % 					graph_matrices.LapMat           : laplacian matrix
    
    n = size(G.Nodes,1);
    graph_edge = G.Edges{:,:};
    adjMat = zeros(n, n);   % Adjacency matrix
    degVec = zeros(n, 1);     % Used for Degree matrix
    incMat = zeros(n, size(graph_edge,1));  % Incidence matrix

    if(isempty(graph_edge) == 0) % graph is not empty
        for ind = 1:size(graph_edge,1) % number of edges
            adjMat(graph_edge(ind,1),graph_edge(ind,2)) = 1;
            adjMat(graph_edge(ind,2),graph_edge(ind,1)) = 1;

            % Each time the vertex appears in edge, increment degVec by 1
            degVec(graph_edge(ind,1),1) = degVec(graph_edge(ind,1),1) + 1;
            degVec(graph_edge(ind,2),1) = degVec(graph_edge(ind,2),1) + 1;

            % Incidence matrix -- rows are the n, and columns are the Edges
            % rowIndex and colIndex: example on K3, 
            % for ind = 1
            %   rowI        = [1, 2], colI = 1; 
            %   incMat(1,1) = 1;    % ij entry corresponding with vertex 1, and edge 12
            %   incMat(2,1) = 1;    % ij entry corresponding with vertex 2, and edge 12
            rowI = graph_edge(ind,:); % row index for ij entry in incMat
            colI = ind; % col index for ij entry in incMat
            incMat(rowI(1), colI) = 1;
            incMat(rowI(2), colI) = 1;
            
        end
        degMat = diag(degVec, 0); % Degree matrix
        LapMat = degMat - adjMat; % Laplacian matrix for undirected graph, assume arbitrary orientation

    else
        degMat = zeros(n);
        LapMat = degMat - adjMat;  
    end

    % Populate data structure graph_matrices
    graph_matrices.adjMat = adjMat;
    graph_matrices.degMat = degMat;
    graph_matrices.incMat = incMat;
    graph_matrices.LapMat = LapMat;
end

