% Genetic algorithm

%PLAN
%{  
    % initialize the population
    % fitness evaluation
    % selection
    % crossover
    % mutation
    % elitism
    % termination and update the population
    % output
%}





function [bestSolution, bestCost, timeTaken] = GeneticAlgorithmCND(NodesFile, LinksFile, LFile, c, runtimeLimit, populationSize, mutationRate, elitismCount)
    % Load the graph data
    Nodes = load(NodesFile);
    Links = load(LinksFile);
    L = load(LFile);
    nNodes = size(Nodes, 1);
    % Create graph
    G = graph(L);
        
    % Initialize population
    population = InitializePopulation(populationSize, nNodes, c);
    timeTaken = 0;
    startTime = tic;
    while toc(startTime) < runtimeLimit       
        % Initialize new population - P'
        newPopulation = zeros(populationSize, c);
        costs = EvaluateCosts(G, population);

        for i = 1:populationSize
            s = Crossover(population, costs, nNodes, c);

            if rand < mutationRate
                s = Mutation(s, nNodes);
            end

            newPopulation(i, :) = s;
        end
        
        % Update population
        population = Selection(G, population, newPopulation, elitismCount);

    end

    costs = EvaluateCosts(G, population);
    [~, bestIdx] = min(costs);
    bestSolution = population(bestIdx, :);
    % cost for the best solution
    bestCost = ConnectedNP(G, bestSolution);
    timeTaken = toc(startTime);
    % Output the best solution
    fprintf("\nGenetic Algorithm for CND Problem\n");
    fprintf("C value: %d\n", c);
    fprintf("Best solution: %s\n", mat2str(bestSolution));
    fprintf("Best cost: %f\n", bestCost);
    
    plotTopology(Nodes, Links, bestSolution);
end

function population = InitializePopulation(populationSize, nNodes, c)
    population = zeros(populationSize, c);
    for i = 1:populationSize
        individual = randperm(nNodes, c);
        while length(unique(individual)) < c
            individual = randperm(nNodes, c);
        end
        population(i, :) = individual;
    end
end

function costs = EvaluateCosts(G, population)
    populationSize = size(population, 1);
    costs = zeros(populationSize, 1);
    for i = 1:populationSize
        costs(i) = ConnectedNP(G, population(i, :));
    end
end

function selected = Selection(G, population, newPopulation, elitismCount)
    % Selection operator:
    % Selects the |P| individuals s from set P U P’ with the best value of f(s)
    % f(s) in this case is the ConnectedNP function (s)
    % Limiting the number of elitist individuals (i.e., elements of P) to a maximum value m (elitismCount)
    
    % all_p = population U newPopulation
    fullPopulation = [population; newPopulation]; 
    % Make sure this is a set by removing duplicates
    fullPopulation = unique(fullPopulation, 'rows');

    % Evaluate the fitness of each individual
    fullCosts = EvaluateCosts(G, fullPopulation);
    
    % Sort the population based on costs
    [~, sortedIdx] = sort(fullCosts);
    
    % Select the best m individuals
    selected = fullPopulation(sortedIdx(1:elitismCount), :);
    
    % Fill the rest of the selected population
    remainingCount = size(population, 1) - elitismCount;
    for i = 1:remainingCount
        competitors = randperm(size(fullPopulation, 1), 2);
        if fullCosts(competitors(1)) < fullCosts(competitors(2))
            selected(elitismCount + i, :) = fullPopulation(competitors(1), :);
        else
            selected(elitismCount + i, :) = fullPopulation(competitors(2), :);
        end
    end
end

% Parent Selection strategies
function parents = TournamentSelection(population, costs)
    % Tournament selection for choosing two parents
    % for each parent, select two random individuals and choose the most fit among the two as the parent
    parents = zeros(2, size(population, 2));
    for i = 1:2
        competitors = randperm(size(population, 1), 2);
        if costs(competitors(1)) < costs(competitors(2))
            parents(i, :) = population(competitors(1), :);
        else
            parents(i, :) = population(competitors(2), :);
        end
    end
end

function parents = FitnessBasedSelection(population, costs)
    % Fitness-based selection for choosing two parents
    % select each parent with a probability proportional to its fitness
    [~, sortedIdx] = sort(costs);
    parents = population(sortedIdx(1:2), :);
end

function parents = RankBasedSelection(population, costs)
    % Rank-based selection for choosing two parents
    % rank all individuals by their fitness value and 
    % randomly select each parent with a probability proportional to its rank

    % Rank the combined population based on costs
    [~, sortedIdx] = sort(costs);

    % Calculate selection probabilities
    n = length(population);
    selectionProbabilities = (1:n) / sum(1:n);

    % Select parents based on probabilities
    parents = zeros(2, size(population, 2));
    for i = 1:2
        r = rand;
        cumulativeProbability = 0;
        for j = 1:n
            cumulativeProbability = cumulativeProbability + selectionProbabilities(j);
            if r <= cumulativeProbability
                selectedIdx = sortedIdx(j);
                break;
            end
        end
        parents(i, :) = population(selectedIdx, :);
    end

    % Return the selected parents
end
% end of Parent Selection strategies

% Crossover operator
function offspring = Crossover(population, costs, nNodes, c)
    % Parent Selection - Choose one of the strategies
    parents = TournamentSelection(population, costs);
    % parents = FitnessBasedSelection(population, costs);
    % parents = RankBasedSelection(population, costs);

    crossoverPoint = randi(c - 1);
    offspring = [parents(1, 1:crossoverPoint), parents(2, crossoverPoint+1:end)];
    offspring = unique(offspring, 'stable');
    
    % Ensure the offspring has exactly c unique nodes
    while length(offspring) < c
        newGene = randi(nNodes);
        if ~ismember(newGene, offspring)
            offspring = [offspring, newGene];
        end
    end
end

function mutatedOffspring = Mutation(offspring, nNodes)
    % One gene of the offspring individual is randomly mutated
    mutationPoint = randi(length(offspring));
    newGene = randi(nNodes);
    while ismember(newGene, offspring)
        newGene = randi(nNodes);
    end
    mutatedOffspring = offspring;
    mutatedOffspring(mutationPoint) = newGene;
end


Nodes_file = 'Nodes200.txt';
Links_file = 'Links200.txt';
L_file = 'L200.txt';

% Genetic algorithm configurations
configs = {
    struct('populationSize', 100, 'mutationRate', 0.5, 'elitismCount', 50),
};


c_values = [8, 10, 12];
runtime_limit = 60; % seconds

% Create results folder if it doesn't exist
results_folder = 'results_time_GA';
if ~exist(results_folder, 'dir')
    mkdir(results_folder);
end

% Loop through each configuration
for k = 1:length(configs)
    config = configs{k};
    populationSize = config.populationSize;
    mutationRate = config.mutationRate;
    elitismCount = config.elitismCount;

    % Run the GA for different values of c and store the results
    for i = 1:length(c_values)
        c = c_values(i);
        results = zeros(10, 2);
        for j = 1:10
            % Run the GA
            [bestSolution, bestCost, elapsedTime] = GeneticAlgorithmCND(Nodes_file, Links_file, L_file, c, runtime_limit, populationSize, mutationRate, elitismCount);
            % Store the best cost and elapsed time in the results matrix
            results(j, 1) = bestCost;
            results(j, 2) = elapsedTime;

            fprintf('Run %d/%d for c=%d completed\n', j, 10, c);
            fprintf("Best cost: %f, Elapsed time: %f\n", bestCost, elapsedTime);

        end

        % Save the results to a file
        resultsFile = sprintf('%s/GRASP_res_%d_%d_%.1f_%d.txt', results_folder, c, populationSize, mutationRate, elitismCount);
        save(resultsFile, 'results', '-ascii', '-tabs');
        fprintf('Results saved to %s\n', resultsFile);
    end
end