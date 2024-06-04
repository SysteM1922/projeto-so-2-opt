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
c = 12;

name = sprintf('CND_lp_%d.lpt', c);
% Open the file for writing
fileID = fopen(name, 'wt');

% Objective function: Minimize the number of connected node pairs
fprintf(fileID, 'min ');
% number of connected node pairs
for i=1:N-1
    for j = i+1 : N
        % u(i,j) = 1 if nodes i and j are connected and 0 otherwise (use the adjacency matrix A)
        fprintf(fileID, ' + u%d_%d', i, j);   
    end
end


fprintf(fileID, '\nsubject to\n');
% constraint 1 -> sum{i=1 to n} 𝑣𝑖 = c 
for i=1:N
    fprintf(fileID, ' + v%d', i);
end
fprintf(fileID, ' = %d\n', c);


% constraint 2 -> 𝑢𝑖𝑗 + 𝑣𝑖 + 𝑣𝑗 ≥ 1 , 𝑖,𝑗 ∈ E

for k=1:A
    i = links(k,1);
    j = links(k,2);
    fprintf(fileID, ' + u%d_%d + v%d + v%d >= 1\n', i, j, i, j);
end

% constraint 3 -> 𝑢𝑖𝑗 ≥ 𝑢{𝑖𝑘} + 𝑢{𝑘𝑗} − 1 + 𝑣𝑘 , 𝑖,𝑗 ∉ 𝐸, 𝑘 ∈ 𝑉(𝑖)
% if i is connected with its neighbour k and k is connected with j, then node i is connected with node j
for i = 1:N
    for j = i+1:N
        
        n_i = find(adj_m(i,:)); % get the neighbours of node i
        
        if any(n_i == j)
            continue;
        end

        % Iterate through all nodes to find common neighbors
        for m = 1:length(n_i)
            k = n_i(m);
            fprintf(fileID, '+ u%d_%d - u%d_%d - u%d_%d - v%d >= -1\n', i, j, min(i, k), max(i, k), min(k, j), max(k, j), k);
        end
    end
end


% Binary constraints
fprintf(fileID, 'binary\n');
for i = 1:N
    fprintf(fileID, 'v%d ', i);
end
fprintf(fileID, '\n');

fprintf(fileID, 'general\n');
for i = 1:N-1
    for j = i+1:N
        fprintf(fileID, 'u%d_%d ', i, j);
    end
end

fprintf(fileID, '\nend');

% Close the file
fclose(fileID);

disp('File generation complete.');