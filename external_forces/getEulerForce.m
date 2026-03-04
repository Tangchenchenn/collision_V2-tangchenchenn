function [Feul, Jeul] = getEulerForce(MultiRod, q, alpha_vec)
% Euler 惯性力：F = -m * (alpha x r)

    n_nodes = MultiRod.n_nodes;
    n_DOF   = MultiRod.n_DOF;

    Feul = zeros(n_DOF, 1);
    Jeul = zeros(n_DOF, n_DOF);

    Alpha = [0, -alpha_vec(3),  alpha_vec(2);
             alpha_vec(3), 0,  -alpha_vec(1);
            -alpha_vec(2), alpha_vec(1), 0];

    for i = 1:n_nodes
        idx = (3*i-2):(3*i);
        m_node = MultiRod.massVec(idx(1));
        pos = q(idx);

        F_node = -m_node * (Alpha * pos);
        Feul(idx) = F_node;

        % 对 q 的导数
        Jeul(idx, idx) = -m_node * Alpha;
    end
end