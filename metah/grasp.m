% Load data
Nodes = load('Nodes200.txt');
Links = load('Links200.txt');
L = load('L200.txt');

% Create Graph
nNodes = size(Nodes, 1);
% disp(nNodes);
G = graph(L);


function [best_solution, best_value, timeTaken] = GRASP(G, c, max_time, r)
    t_start = tic;

    solution = GreedyRand(G, c, r);
    [best_value, best_solution, duration] = steepestAscentHillClimbing(G, solution, @bestNeighbor1);
    timeTaken = duration;
    while toc(t_start) < max_time
        % shuffle seed
        rng('shuffle'); 

        % Constructive Phase (Greedy Randomized)
        solution = GreedyRand(G, c, r);
        
        % Search Phase (Adaptive Search - SA-HC)
        [value, solution, duration] = steepestAscentHillClimbing(G, solution, @bestNeighbor1);
        % [value, solution, duration] = steepestAscentHillClimbing(G, solution, @bestNeighbor2);
        %fprintf("Solution: %s, Value: %f\n", mat2str(solution), value);
        % Evaluate the solution
        if value < best_value
            best_solution = solution;
            best_value = value;
            timeTaken = duration;
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
        e = R(randi(min(r, size(R, 1))), 1); % Ensure r is within bounds
        sol = [sol e];
        E = setdiff(E, e);
    end

    %fprintf('Greedy Randomized solution: %s\n', mat2str(sol));
end

% SA-HC function
function [nodesConnected, bestSelected, duration] = steepestAscentHillClimbing(G, bestSelected, selectFunction)
    nodesConnected = inf;
    improved = true;
    t = tic;
    while improved
        improved = false;
        [bestNeighbor, connected] = selectFunction(G, bestSelected);
        if connected < nodesConnected
            nodesConnected = connected;
            bestSelected = bestNeighbor;
            improved = true;
            % fprintf('Improved in local search: %s with cost %f\n', mat2str(bestSelected), nodesConnected);
        end

    end
    % fprintf('Best solution: %s, cost: %f\n', mat2str(solution), nodesConnected);

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



%{
    % different settings
    c_values = [8, 10, 12];
    max_time = 60; % seconds
    r_values = [3, 5, 10, 50]; % Different r values for GreedyRand

    for r = r_values
        results = zeros(10, length(c_values)); % Initialize results matrix for 10 runs and 3 c values
        for c_idx = 1:length(c_values)
            c = c_values(c_idx);
            for run = 1:10
                [~, best_value, timeTaken] = GRASP(G, c, max_time, r);
                results(run, c_idx) = best_value;
            end
        end
        
        % Write results to file
        filename = sprintf('GRASP2_RES_%d.txt', r);
        fileID = fopen(filename, 'w');
        fprintf(fileID, 'Results for r = %d\n', r);
        fprintf(fileID, 'C8\tC10\tC12\n');
        for run = 1:10
            fprintf(fileID, '%f\t%f\t%f\n', results(run, 1), results(run, 2), results(run, 3));
        end
        fclose(fileID);
    end 
%}

% best settings 
r = 50;
c_values = [8, 10, 12];
max_time = 60; % seconds
% Create results folder if it doesn't exist
results_folder = 'results_time_GRASP';
if ~exist(results_folder, 'dir')
    mkdir(results_folder);
end

for c = c_values
    results = zeros(10, 2); % Initialize results matrix for 10 runs
    for run = 1:10
        [~, best_value, timeTaken] = GRASP(G, c, max_time, r);
        results(run, 1) = best_value;
        results(run, 2) = timeTaken;

        fprintf('Run %d/%d for c=%d completed\n', run, 10, c);
        fprintf("Best value: %f, Time taken: %f\n", best_value, timeTaken)

    end

    % Write results to file
    filename = sprintf('%s/GRASP_res_%d_%d.txt', results_folder, c, r);
    save(filename, 'results', '-ascii', '-tabs');
    fprintf('Results saved to %s\n', filename);
end


