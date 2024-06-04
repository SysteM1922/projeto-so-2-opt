function out=AverageSP(G,selected)
% AverageSP(G,servers) - Computes the average shortest path length from each
%         node to its closest selected node (returns -1 for invalid input data)
%
% G:         graph of the network
% selected:  a row array with IDs of selected nodes
    
    if length(selected)<1
        out= -1;
        return
    end
    nNodes= numnodes(G);
    if (max(selected)>nNodes || min(selected)<1 || length(unique(selected))<length(selected))
        out= -1;
        return
    end
    aux= setdiff(1:nNodes,selected);
    dist= distances(G,selected,aux);
    if length(selected)>1
        out= sum(min(dist))/nNodes;
    else
        out= sum(dist)/nNodes;
    end
end