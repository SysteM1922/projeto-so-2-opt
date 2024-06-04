% Load data
Nodes = load('Nodes200.txt');
Links = load('Links200.txt');
L = load('L200.txt');

% Create Graph
nNodes = size(Nodes, 1);
G = graph(L);

% GRASP function
function [best_solution, best_value] = GRASP(G, c, max_time)
    best_solution = [];
    best_value = inf;
    t_start = tic;

    while toc(t_start) < max_time
        % Greedy Randomized Construction Phase
        solution = GreedyRandomized(G, c);
        
        % Adaptive Search Phase using VND
        [solution, value] = VND(G, solution);

        if value < best_value
            best_solution = solution;
            best_value = value;
        end
    end
end

% Greedy Randomized function
function solution = GreedyRandomized(G, c)
    nNodes = numnodes(G);
    solution = [];
    
    % Greedy construction
    while length(solution) < c
        candidates = setdiff(1:nNodes, solution); % Nodes not yet in the solution
        candidate_values = arrayfun(@(x) ConnectedNP(G, [solution x]), candidates); % Evaluate candidates using ConnectedNP
        [~, best_idx] = min(candidate_values); % Find the best candidate
        
        % Select the best candidate
        new_node = candidates(best_idx);
        
        solution = [solution new_node]; % Add the chosen node to the solution
    end
end

% Variable Neighborhood Descent (VND) function
function [solution, best_value] = VND(G, solution)
    kmax = 2; % Number of neighborhood structures
    k = 1;
    best_value = ConnectedNP(G, solution);

    while k <= kmax
        if k == 1
            % Neighborhood structure 1: Swap nodes
            [new_solution, new_value] = SwapNodes(G, solution, best_value);
        elseif k == 2
            % Neighborhood structure 2: Replace least connected node
            [new_solution, new_value] = ReplaceLeastConnected(G, solution, best_value);
        end
        
        if new_value < best_value
            solution = new_solution;
            best_value = new_value;
            k = 1;
        else
            k = k + 1;
        end
    end
end

% Neighborhood structure 1: Swap nodes
function [best_neighbor, best_value] = SwapNodes(G, solution, current_value)
    best_neighbor = solution;
    best_value = current_value;
    nNodes = numnodes(G);

    for i = 1:length(solution)
        for j = setdiff(1:nNodes, solution)
            neighbor = solution;
            neighbor(i) = j;

            value = ConnectedNP(G, neighbor);
            if value < best_value
                best_value = value;
                best_neighbor = neighbor;
            end
        end
    end
end

% Neighborhood structure 2: Replace least connected node
function [best_neighbor, best_value] = ReplaceLeastConnected(G, solution, current_value)
    best_neighbor = solution;
    best_value = current_value;
    nNodes = numnodes(G);

    % Find the node with the least neighbors in the current solution
    min_neighbors = inf;
    node_to_swap = -1;
    for i = 1:length(solution)
        num_neighbors = degree(G, solution(i));
        if num_neighbors < min_neighbors
            min_neighbors = num_neighbors;
            node_to_swap = i;
        end
    end
    
    % Swap with a node outside the solution
    candidates = setdiff(1:nNodes, solution);
    for j = candidates
        neighbor = solution;
        neighbor(node_to_swap) = j;

        value = ConnectedNP(G, neighbor);
        if value < best_value
            best_value = value;
            best_neighbor = neighbor;
        end
    end
end

% ConnectedNP function (for reference)
function out = ConnectedNP(G, selected)
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

% Run GRASP with VND
c = 12;
max_time = 60; % seconds

[best_solution, best_value] = GRASP(G, c, max_time);
fprintf('C: %d\n', c);
fprintf('Best solution: %s, cost: %f\n', mat2str(best_solution), best_value);
fprintf('---------------------------------------\n');
