% Load data
Nodes = load('Nodes200.txt');
Links = load('Links200.txt');
L = load('L200.txt');

% Create Graph
nNodes = size(Nodes, 1);
% disp(nNodes);
G = graph(L);


function [best_solution, best_value] = GRASP(G, c, max_time, r, local_search_time)
    best_solution = [];
    best_value = inf;
    t_start = tic;

    while toc(t_start) < max_time
        % Greedy Randomized Construction Phase
        % solution = GreedyRandomized(G, c);
        
        solution = GreedyRand(G, c, r);

        % Adaptive Search Phase (SA-HC)
        [solution, value] = AdaptiveSearch(G, solution);

        % Adaptive Search Phase (SA-HC) with heuristic
        % [solution, value] = AdaptiveSearchHeuristic(G, solution);
        
        % SA-HC function GUI
        % [solution, value, ~, ~] = steepestAscentHillClimbing(G, solution, @bestNeighbor1, local_search_time);

        % Evaluate the solution
        if value < best_value
            best_solution = solution;
            best_value = value;
        end
    end
end

% Greedy Randomized function - not actually randomized
%{
  pick an empty solution
  for i = 1 to c
    pick a set of nodes that looks best
    select 1 random from it and add it to the solution
    remove the selected node from the set
%}
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

    % fprintf('Greedy Randomized solution: %s\n', mat2str(solution));
end

% from the slides
function sol = GreedyRand(G, c, r)
    E = 1:numnodes(G);
    sol = [];
    for i = 1:c
        R = [];
        for j=E
            R = [ R ; j ConnectedNP(G, [sol j]) ];
        end
        R = sortrows(R, 2);
        e = R(randi(r), 1);
        sol = [sol e];
        E = setdiff(E, e);
    end

    % fprintf('Greedy Randomized solution: %s\n', mat2str(sol));
end

% SA-HC function GUI with time limit
function [nodesConnected, bestSelected, iteractions, duration] = steepestAscentHillClimbing(G, bestSelected, selectFunction, max_time)
    nodesConnected = inf;
    improved = true;
    iteractions = 0;
    t = tic;

    while improved && toc(t) < max_time
        improved = false;
        iteractions = iteractions + 1;
        [bestNeighbor, connected] = selectFunction(G, bestSelected);
        if connected < nodesConnected
            nodesConnected = connected;
            bestSelected = bestNeighbor;
            improved = true;
        end
    end

    duration = toc(t);
end

% a neighbor solution is obtained by swapping a node in the current solution by a node not in the current solution
function [bestNeighbor, bestConnected] = bestNeighbor1(G, current)
    nNodes = numnodes(G);
    others = setdiff(1:nNodes, current);
    bestConnected = inf;
    bestNeighbor = current;
    for a = current
        for b = others
            neighbor = [setdiff(current, a), b];
            connected = ConnectedNP(G, neighbor);
            if connected < bestConnected
                bestConnected = connected;
                bestNeighbor = neighbor;
            end
        end
    end
end

% Adaptive Search (SA-HC) function
function [solution, best_neighbor_value] = AdaptiveSearch(G, solution)
    improved = true;
    while improved
        improved = false;
        nNodes = numnodes(G);
        best_neighbor_value = ConnectedNP(G, solution);
        best_neighbor = solution;

        for i = 1:length(solution)
            for j = setdiff(1:nNodes, solution)
                neighbor = solution;
                neighbor(i) = j;
                
                value = ConnectedNP(G, neighbor);
                if value < best_neighbor_value
                    best_neighbor_value = value;
                    best_neighbor = neighbor;
                    improved = true;
                    % fprintf('Improved in local search: %s with cost %f\n', mat2str(best_neighbor), best_neighbor_value);
                end
            end
        end

        if improved
            solution = best_neighbor;
        end
    end
end

% Adaptive Search function (with heuristic of least neighbors)
% choose the node with the least neighbors in the current solution to swap
% for some reason, not working the best (maybe try to swap it with a neighbor node instead of any node?)
% increase time lmit just to see if it works
function [solution, best_neighbor_value] = AdaptiveSearchHeuristic(G, solution)
    improved = true;
    while improved
        improved = false;
        best_neighbor_value = ConnectedNP(G, solution); 
        best_neighbor = solution;

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
        nNodes = numnodes(G);
        candidates = setdiff(1:nNodes, solution);
        for j = candidates
            neighbor = solution;
            neighbor(node_to_swap) = j;

            value = ConnectedNP(G, neighbor); 
            if value < best_neighbor_value
                best_neighbor_value = value;
                best_neighbor = neighbor;
                improved = true;
                % fprintf('Improved in local search: %s with cost %f\n', mat2str(best_neighbor), best_neighbor_value);
            end
        end

        if improved
            solution = best_neighbor;
        end
    end
end

% Experiment with different settings and both adaptive search functions
c_values = [8, 10, 12];
max_time = 60; % seconds

r = 3;  % r value for GreedyRand
local_search_time = 13; % seconds for local search (SA-HC) - 13 s is just an example, still need to find the best value!


for c = c_values
    [best_solution, best_value] = GRASP(G, c, max_time, r, local_search_time);
    fprintf('C: %d\n', c);
    fprintf('Best solution: %s, cost: %f\n', mat2str(best_solution), best_value);
    fprintf('---------------------------------------\n');
end

%{
% run for different c values 10 times and print the results (min, avg, max)
results = struct('c', {}, 'min_value', {}, 'avg_value', {}, 'max_value', {});

for i = 1:length(c_values)
    
    c = c_values(i);
    values = zeros(1, 10);
    
    for k = 1:10
        [~, values(k)] = GRASP(G, c, max_time);
    end
    
    results(end+1).c = c;
    results(end+1).min_value = min(values);
    results(end+1).avg_value = mean(values);
    results(end+1).max_value = max(values);
end

% Print results
for i = 1:length(results)
    fprintf('C: %d\n', results(i).c);
    fprintf('Min Value: %f, Avg Value: %f, Max Value: %f\n', results(i).min_value, results(i).avg_value, results(i).max_value);
    fprintf('---------------------------------------\n');
end

%}