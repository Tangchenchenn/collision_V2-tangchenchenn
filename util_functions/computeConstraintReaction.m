function reaction = computeConstraintReaction(MultiRod, q, f, sim_params)
% 从 fixedDOF 的残差提取约束反力/反力矩
%
% 这里的 f 是 timeStepper 收敛时的“全量残差”
% 对约束自由度，有：
%   M a - F_internal - F_external - R_constraint = 0
% 因此：
%   R_constraint = f
%
% motor 需要提供的扭矩 = - 约束对绳子的反力矩

    reaction = struct();

    fixedDOF = MultiRod.fixedDOF(:);
    reaction.full = zeros(size(f));
    reaction.full(fixedDOF) = f(fixedDOF);

    if isfield(sim_params, 'hub_center')
        hub_center = sim_params.hub_center(:);
    else
        hub_center = [0; 0; 0];
    end

    fixed_nodes = MultiRod.fixed_nodes(:)';
    n_fix = numel(fixed_nodes);

    reaction.node_force   = zeros(3, n_fix);
    reaction.node_moment  = zeros(3, n_fix);

    for k = 1:n_fix
        node_id = fixed_nodes(k);
        idx = mapNodetoDOF(node_id);

        Fi = reaction.full(idx);          % 约束对绳子的反力
        ri = q(idx) - hub_center;         % 该固定点相对轮毂中心的位置

        reaction.node_force(:, k)  = Fi;
        reaction.node_moment(:, k) = cross(ri, Fi);
    end

    reaction.total_force  = sum(reaction.node_force, 2);
    reaction.total_moment = sum(reaction.node_moment, 2);

    % 对绳子的约束反力矩
    reaction.constraint_torque_z = reaction.total_moment(3);

    % 电机为维持给定 omega，需要输出相反的扭矩
    reaction.motor_torque_z = -reaction.constraint_torque_z;

    % 预留：若以后固定 edge theta DOF，可继续补 fixed edge 的扭转反力
    reaction.fixed_edge_torque = [];
end