function x_out = element_graph_dilation(x_in, A)
% ELEMENT_GRAPH_DILATION  Graph dilation: each node takes maximum over neighborhood+self
N = numel(x_in);
x_out = x_in;
for i = 1:N
    nb = find(A(i,:));
    if isempty(nb), idx = i; else, idx = [i, nb]; end
    x_out(i) = max(x_in(idx));
end
end
