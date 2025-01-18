function [W_opt, fval, err] = opadmm_gl(X, alpha, beta, W_ref, param_solver)

    % parse input parameters
    gamma         = param_solver.gamma;
    rho_0         = param_solver.rho_0;
    tau1_fun      = param_solver.tau1_fun;
    tau2_fun      = param_solver.tau2_fun;
    inner_loops   = param_solver.inner_loops;
    th            = param_solver.th;

    % initialization
    m = size(X,1);
    n = size(X,2);
    p = m*(m-1)/2;

    idx = 1;
    n_blocks = length(W_ref);
    n_samp_per_block = n/n_blocks;

    fval = zeros(n,1);
    err  = zeros(n,1);

    % S initialization
    [S,St] = sum_squareform(m);

    % variables initilization
    z_k   = zeros(p,1);
    w_k   = zeros(p,1);
    Sw    = zeros(m,1);
    v_k   = Sw;
    lmb_k = zeros(m,1);

    W_opt = zeros(p,n);

    rho_k  = rho_0;
    tau1_k = tau1_fun(rho_k);
    tau2_k = tau2_fun(rho_k);

    gamma_k = gamma;

    t_0 = tic;

    for k = 1:n

        if k == (idx-1)*n_samp_per_block+1
            w_ref = W_ref{idx};
            idx = idx+1;
        end

        if n_blocks == 1
            gamma_k = 1/k;
        end

        % online adaption
        z_km1 = z_k;
        Z     = sparse(gsp_distanz(X(:,k)').^2);
        z     = squareform_sp(Z);
        z_k   = gamma_k*z + (1-gamma_k)*z_km1;

        for r = 1:inner_loops

            % primary variable update
            w_km1 = w_k;
            w_l = w_k - tau1_k*(rho_k*St*(Sw-v_k) + St*lmb_k + 2*z_k);
            w_k = max(w_l,0)/(2*tau1_k*beta+1);

            Sw = S*w_k;

            % auxiliary variable update
            v_km1 = v_k;
            v_l   = (1-tau2_k*rho_k)*v_k + tau2_k*rho_k*Sw + tau2_k*lmb_k;
            v_k   = (v_l + sqrt(v_l.^2+4*tau2_k*alpha))/2;

            % dual variable update
            lmb_km1 = lmb_k;
            d_k     = Sw - v_k;
            lmb_k   = lmb_km1 + rho_k*d_k;

        end

        fprintf('Iteration %d completed! Elapsed time: %.2f sec.\n', k, toc(t_0));

        logv_k = log(v_k);
        logv_k(v_k == 0) = 0;
 
        fval(k) = 2*z_k'*w_k - alpha*sum(logv_k) + beta*(w_k'*w_k);

        w_opt_k               = w_k;
        w_opt_k(w_opt_k < th) = 0;

        W_opt(:,k) = w_opt_k;

        err(k) = norm(w_opt_k-w_ref);

    end

end