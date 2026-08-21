function mdl = assign_electrode_nodes(mdl, node_ids)
%ASSIGN_ELECTRODE_NODES  Assign node indices to electrodes in a forward model
%   mdl = assign_electrode_nodes(mdl, node_ids)
%   node_ids: vector of node indices (one per electrode)

for idx = 1:numel(node_ids)
    mdl.electrode(idx).nodes = node_ids(idx);
end
end
