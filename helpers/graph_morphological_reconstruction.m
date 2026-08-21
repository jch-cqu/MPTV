function recon = graph_morphological_reconstruction(marker, mask, A)
% GRAPH_MORPHOLOGICAL_RECONSTRUCTION  Graph-based geodesic reconstruction
% recon = graph_morphological_reconstruction(marker, mask, A)
% marker, mask are logical vectors (Ne x 1); A is sparse adjacency (Ne x Ne)
N = numel(marker);
cur = logical(marker(:));
prev = false(size(cur));
while ~isequal(cur, prev)
    prev = cur;
    expand = false(size(cur));
    for i = 1:N
        if mask(i)
            nb = find(A(i,:));
            if any(cur(nb)) || cur(i)
                expand(i) = true;
            end
        end
    end
    cur = cur | expand;
end
recon = cur;
end
