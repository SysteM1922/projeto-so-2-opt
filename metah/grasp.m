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
        % shuffle seed
        rng('shuffle'); 

        % Constructive Phase (Greedy Randomized)
        solution = GreedyRand(G, c, r);
        
        % Search Phase (Adaptive Search - SA-HC)
        [value, solution, ~, ~] = steepestAscentHillClimbing(G, solution, @bestNeighbor1, local_search_time);
        % [value, solution, ~, ~] = steepestAscentHillClimbing(G, solution, @bestNeighbor2, local_search_time);
        %fprintf("Solution: %s, Value: %f\n", mat2str(solution), value);
        % Evaluate the solution
        if value < best_value
            best_solution = solution;
            best_value = value;
        end
    end
end

% Greedy Randomized function 
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
        % e = R(randi(min(r, size(R, 1))), 1); % Ensure r is within bounds
        sol = [sol e];
        E = setdiff(E, e);
    end

    %fprintf('Greedy Randomized solution: %s\n', mat2str(sol));
end

% SA-HC function GUI with time limit 
function [nodesConnected, bestSelected, iteractions, duration] = steepestAscentHillClimbing(G, bestSelected, selectFunction, max_time)
    nodesConnected = inf;
    improved = true;
    iteractions = 0;
    t = tic;
    solution = bestSelected;

    while improved && toc(t) < max_time
        improved = false;
        iteractions = iteractions + 1;
        [bestNeighbor, connected] = selectFunction(G, bestSelected);
        if connected < nodesConnected
            nodesConnected = connected;
            bestSelected = bestNeighbor;
            improved = true;
            % fprintf('Improved in local search: %s with cost %f\n', mat2str(bestSelected), nodesConnected);
        end

    end

    solution = bestSelected;
    %fprintf('Best solution: %s, cost: %f\n', mat2str(solution), nodesConnected);

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

% a neighbor solution is obtained by swapping a node in the current solution by a neighbor node not in the current solution
function [bestNeighbor, bestConnected] = bestNeighbor2(G,current)
    bestConnected = inf;
    bestNeighbor = current;
    for a = current
        others = setdiff(neighbors(G,a)',current);  % get the neighbors of the current node that are not in the current solution
        for b = others
            neighbor = [setdiff(current, a), b];    % swap a node in the current solution by a neighbor node not in the current solution
            connected = ConnectedNP(G,neighbor);
            if connected < bestConnected
                bestConnected = connected;
                bestNeighbor = neighbor;
            end
        end
    end
end

% ConnectedNP function
function out=ConnectedNP(G,selected)
    % ConnectedNP(G,selected) - Computes the number of node pairs that can communicate
    %         if the selected nodes are eliminated (returns -1 for invalid input data)
    %
    % G:         graph of the network
    % selected:  a row array with IDs of selected nodes
        
    nNodes= numnodes(G);
    if length(selected)>=1
        if (max(selected)>nNodes || min(selected)<1 || length(unique(selected))<length(selected))
            out= -1;
            return
        end
    end
    aux= setdiff(1:nNodes,selected);
    Gr= subgraph(G,aux);
    dist= distances(Gr);
    out= (sum(dist(:)<Inf) - numnodes(Gr))/2;
end

% Experiment with different settings and both adaptive search functions
c_values = [8, 10, 12];
max_time = 60; % seconds

r = 3;  % r value for GreedyRand
local_search_time = 13; % seconds for local search (SA-HC) - 13 s is just an example, still need to find the best value!



%{
 for c = c_values
    [best_solution, best_value] = GRASP(G, c, max_time, r, local_search_time);
    fprintf('C: %d\n', c);
    fprintf('Best solution: %s, cost: %f\n', mat2str(best_solution), best_value);
    fprintf('---------------------------------------\n');
    end
%}
%{
    c = 8;
    [best_solution, best_value] = GRASP(G, c, max_time, r, local_search_time);
    fprintf('C: %d\n', c);
    fprintf('Best solution: %s, cost: %f\n', mat2str(best_solution), best_value);
    fprintf('---------------------------------------\n');
%}

% run for different c values 10 times and print the results (min, avg, max)
results = struct('c', {}, 'min_value', {}, 'avg_value', {}, 'max_value', {});

for i = 1:length(c_values)
    
    c = c_values(i);
    values = zeros(1, 10);
    
    for k = 1:10
        [best_solution, best_value] = GRASP(G, c, max_time, r, local_search_time);
        values(k) = best_value;
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

