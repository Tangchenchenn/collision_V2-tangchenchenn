function [environment,imc] = createEnvironmentAndIMCStructs(env,geom,material,sim_params)
environment = struct();
imc = struct();

environment.ext_force_list = env.ext_force_list;

if ismember("gravity",env.ext_force_list)
    environment.g = env.g;
end

if ismember("buoyancy",env.ext_force_list)
    environment.rho = env.rho;
end

if ismember("viscous", env.ext_force_list)
    environment.eta = env.eta;
end

if ismember("aerodynamic", env.ext_force_list)
    environment.rho = env.rho;
    environment.Cd  = env.Cd;
    environment.rod_r0 = geom.rod_r0;   % <<< 新增
end

if ismember("pointForce", env.ext_force_list)
    environment.ptForce = env.ptForce;
    environment.ptForce_node = env.ptForce_node;
end

if ismember("rft", env.ext_force_list)
    environment.ct = env.ct;
    environment.cn = env.cn;
end

if ismember("selfContact", env.ext_force_list)
    imc.k_c = material.contact_stiffness;
    imc.contact_len = 2 * geom.rod_r0;
    imc.delta = 0.5 * 0.3125 * imc.contact_len;
    imc.omega = 20;
    imc.scale = 1 / geom.rod_r0;
    imc.C = [];

    if isfield(env, 'contact_params')
        imc.ice_radius = env.contact_params.ice_radius;

        if isfield(env.contact_params, 'num_ice')
            imc.num_ice = env.contact_params.num_ice;
            imc.array_radius = env.contact_params.array_radius;
            imc.array_center_dist = env.contact_params.array_center_dist;
            imc.omega_spin = env.contact_params.omega_spin;
            imc.peak_force = env.contact_params.peak_force;
        else
            imc.ice_center_dist = env.contact_params.ice_center_dist;
        end

        imc.rod_radius = env.contact_params.rod_radius;
        imc.omega_mag  = env.contact_params.omega_mag;
        imc.sigma_t    = env.contact_params.sigma_t;
        imc.z_root     = env.contact_params.z_root;
        imc.is_broken  = env.contact_params.is_broken;
        imc.ice_z_min  = env.contact_params.ice_z_min;
        imc.ice_z_max  = env.contact_params.ice_z_max;
        imc.theta_accumulated = 0;

        if isfield(env.contact_params, 'active_time')
            imc.active_time = env.contact_params.active_time;   % <<< 新增
        end
    end

    % <<< 修改：既支持 selfFriction，也支持 contact_params.compute_friction
    want_friction = ismember("selfFriction", env.ext_force_list);

    if isfield(env, 'contact_params') && isfield(env.contact_params, 'compute_friction')
        want_friction = want_friction || logical(env.contact_params.compute_friction);
    end

    if want_friction
        imc.compute_friction = true;
        imc.mu_k = material.mu;
        imc.velTol = material.velTol;
    else
        imc.compute_friction = false;
        imc.mu_k = 0;
        imc.velTol = 0;
    end
end

if ismember("floorContact", env.ext_force_list)
    imc.floor_z = env.floor_z;
    imc.k_c_floor = env.contact_stiffness;
    imc.h = geom.rod_r0;
    imc.delta_floor = 1 * geom.rod_r0;
    imc.omega = 20;
    imc.scale = 1 / imc.h;
    environment.showFloor = true;

    if ismember("floorFriction", env.ext_force_list)
        imc.floor_has_friction = true;
        imc.mu_floor = env.mu;
        imc.velTol = env.velTol;
    else
        imc.floor_has_friction = false;
    end
end

if sim_params.static_sim
    environment.static_g = env.g;
end
end