function x_out = element_graph_erosion(x_in, A)
% ELEMENT_GRAPH_EROSION  Graph erosion: each node takes minimum over neighborhood+self
N = numel(x_in);
x_out = x_in;
for i = 1:N
    nb = find(A(i,:));
    if isempty(nb), idx = i; else, idx = [i, nb]; end
    x_out(i) = min(x_in(idx));
end
end
