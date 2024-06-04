% exact method based on integer linear programming
nodes = load('Nodes200.txt');
links = load('Links200.txt');
lengths = load('L200.txt');
L = load('L200.txt');

% Create the graph
G = graph(L);

% get the adjacency matrix
adj_m = adjacency(G);

N = 200;    % number of nodes (nodes)
A = 250;    % number of links (L)

% c = 8, 10 and 12 (number of critical nodes)
settings = [12, 10, 8];  % Different values of c

% NOTES:
%{
  
    % GRASP -> Greedy Randomized Adaptive Search Procedure 
    % two phases: construction and local search.
    % The construction phase builds a feasible solution, and the local search phase attempts to improve it.
    % Each iteration provides a candidate solution, and the best one is selected.

    % 1 Generate initial random solution (S) with c critical nodes
    % 2 Construction phase: 
        % Create a candidate solution by iteratively adding nodes based on a greedy function.
        % The greedy function selects the node that maximizes the number of connected node pairs.

    % 3 Local search phase:
        % Improve the candidate solution by exploring the neighborhood (e.g., swapping nodes).
        % Improve the candidate solution by iteratively removing and adding nodes.
        % The local search phase uses a greedy function to select the node to remove or add.

    % Stop criterion: 
        % Set by runtime limit (60 seconds per run).
%}

% GRASP parameters
max_runtime = 60; % runtime limit in seconds

for c = settings
    fprintf('Testing with c = %d and max_runtime = %d\n', c, max_runtime);
    % Initialize the timer and best solution for this setting
    start_time = tic;
    best_solution = [];
    best_cost = Inf;
    
    % Main loop
    while toc(start_time) < max_runtime
        % Greedy Randomized: Generate initial random solution with c critical nodes
        S = greedy_randomized(G, c, N);
    
        % Adaptive Search: Improve the candidate solution
        % [S, cost] = adaptive_search(G, S, c);

        % SAHC: Improve the candidate solution using Steepest Ascent Hill Climbing
        [S, cost] = SAHC(G, S, c);
        
        % VND: Improve the candidate solution using VND
        % [S, cost] = VND(G, S, c);

        % Update best solution if a better one is found
        if cost < best_cost
            best_solution = S;
            best_cost = cost;
        end
    end

    % Display the best solution
    fprintf('Best solution found:\n');
    disp(best_solution);
    fprintf('Cost of best solution: %f\n', best_cost);

end    


function S = greedy_randomized(G, c, N)
    % Generate an initial random solution
    available_nodes = 1:N;
    random_indices = randperm(N, c);
    S = available_nodes(random_indices);

    % Print the initial random solution
    %fprintf('Initial random solution: %s\n', mat2str(S));

    % Improve the solution by a greedy function
    for i = 1:c
        best_node = -1;
        best_value = -inf;
        for node = available_nodes
            if ~ismember(node, S)
                temp_solution = [S node];
                % value = greedy_function(G, temp_solution);
                value = greedy_function1(G, temp_solution);
                if value > best_value
                    best_value = value;
                    best_node = node;
                end
            end
        end
        if best_node ~= -1
            S(i) = best_node;
        end
    end

    % Print the improved solution
    %fprintf('Improved solution: %s\n', mat2str(S));
end

function value = greedy_function(G, S)
    % The greedy function selects the node that maximizes the number of connected node pairs.
    subG = subgraph(G, S);
    value = numedges(subG);
end

function value = greedy_function1(G, S)
    subG = subgraph(G, S);
    components = conncomp(subG);
    largest_component_size = max(histcounts(components));
    value = -largest_component_size;  % Prioritize smaller largest component
end

function [S, cost] = VND(G, S, c)
    % Improve the candidate solution using VND
    kmax = 2; % Number of neighborhood structures
    k = 1;
    best_cost = ConnectedNP(G, S);
    s_best = S;

    while k <= kmax
        s = local_search(G, s_best, c, k);
        s_cost = ConnectedNP(G, s);

        if s_cost < best_cost
            s_best = s;
            best_cost = s_cost;
            k = 1; % Return to the first neighborhood structure
        else
            k = k + 1; % Move to the next neighborhood structure
        end
    end

    S = s_best;
    cost = best_cost;
end

function s = local_search(G, S, c, k)
    % Local search for the given neighborhood structure k
    % Neighborhood structures:
    % k = 1: Swap one node
    % k = 2: Swap two nodes

    best_cost = ConnectedNP(G, S);
    s = S;

    if k == 1
        % Swap one node
        for i = 1:c
            for node = 1:numnodes(G)
                if ~ismember(node, S)
                    new_solution = S;
                    new_solution(i) = node;
                    new_cost = ConnectedNP(G, new_solution);
                    if new_cost < best_cost
                        s = new_solution;
                        best_cost = new_cost;

                        % Print the improved solution during local search
                        fprintf('{1} Improved in local search: %s with cost %f\n', mat2str(s), best_cost);
                    end
                end
            end
        end
    elseif k == 2
        % Swap two nodes
        for i = 1:c
            for j = i+1:c
                for node1 = 1:numnodes(G)
                    for node2 = node1+1:numnodes(G)
                        if ~ismember(node1, S) && ~ismember(node2, S)
                            new_solution = S;
                            new_solution([i, j]) = [node1, node2];
                            new_cost = ConnectedNP(G, new_solution);
                            if new_cost < best_cost
                                s = new_solution;
                                best_cost = new_cost;

                                % Print the improved solution during local search
                                fprintf('{2} Improved in local search: %s with cost %f\n', mat2str(s), best_cost);
                            end
                        end
                    end
                end
            end
        end
    end
end

function [S, cost] = adaptive_search(G, S, c)
    % Improve the candidate solution by exploring the neighborhood (e.g., swapping nodes).
    best_cost = ConnectedNP(G, S);
    improvement = true;

    while improvement
        improvement = false;
        % Try removing and adding nodes
        for i = 1:c
            for node = 1:numnodes(G)
                if ~ismember(node, S)
                    new_solution = S;
                    new_solution(i) = node;
                    new_cost = ConnectedNP(G, new_solution);
                    if new_cost < best_cost
                        S = new_solution;
                        best_cost = new_cost;
                        improvement = true;

                        % Print the improved solution during local search
                        fprintf('Improved in local search: %s with cost %f\n', mat2str(S), best_cost);
                    end
                end
            end
        end
    end

    cost = best_cost;
end

function [S, cost] = SAHC(G, S, c)
    % Improve the candidate solution using Steepest Ascent Hill Climbing
    improvement = true;
    best_cost = ConnectedNP(G, S);

    while improvement
        improvement = false;
        best_neighbor = S;
        best_neighbor_cost = best_cost;

        % Explore the neighborhood by swapping each node in S with a node not in S
        for i = 1:c
            for node = 1:numnodes(G)
                if ~ismember(node, S)
                    new_solution = S;
                    new_solution(i) = node;
                    new_cost = ConnectedNP(G, new_solution);
                    if new_cost < best_neighbor_cost
                        best_neighbor = new_solution;
                        best_neighbor_cost = new_cost;
                        improvement = true;
                        % fprintf('Improved in local search: %s with cost %f\n', mat2str(best_neighbor), best_neighbor_cost);
                        % Break the inner loop to explore the next neighbor
                    end
                end
            end
        end

        % Update the current solution to the best neighbor
        if improvement
            S = best_neighbor;
            best_cost = best_neighbor_cost;
        end
    end

    cost = best_cost;
end


function out = ConnectedNP(G, selected)
    % ConnectedNP(G,selected) - Computes the number of node pairs that can communicate
    % if the selected nodes are eliminated (returns -1 for invalid input data)
    %
    % G: graph of the network
    % selected: a row array with IDs of selected nodes
    
    nNodes = numnodes(G);
    if length(selected) >= 1
        if (max(selected) > nNodes || min(selected) < 1 || length(unique(selected)) < length(selected))
            out = -1;
            return;
        end
    end
    aux = setdiff(1:nNodes, selected);
    Gr = subgraph(G, aux);
    dist = distances(Gr);
    out = (sum(dist(:) < Inf) - numnodes(Gr)) / 2;
end