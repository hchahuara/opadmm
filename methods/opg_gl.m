function [W_opt, fval, err] = opg_gl(X, alpha, beta, W_ref, param_solver)

    % parse input parameters
    gamma       = param_solver.gamma;
    inner_loops = param_solver.inner_loops;
    th          = param_solver.th;

    % initialization
    m = size(X,1);
    n = size(X,2);
    p = m*(m-1)/2;

    idx = 1;
    n_blocks = length(W_ref);
    n_samp_per_block = n/n_blocks;

    z_k = zeros(p,1);
    w_k = zeros(p,1);

    W_opt = zeros(p,n);

    fval = zeros(n,1);
    err  = zeros(n,1);

    [S,St] = sum_squareform(m);

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

            Sw = S*w_k;

            grad_k = 2*beta*w_k - alpha*St*(1./(Sw+eps));
            mu_k   = 1/(2*beta + 2*alpha*(m-1)/(min(Sw)^2));
 
            w_l = w_k - mu_k*(grad_k+2*z_k);

            w_k = max(eps, w_l);
            w_t = max(0, w_l);

        end

        fprintf('Iteration %d completed! Elapsed time: %.2f sec.\n', k, toc(t_0));

        Sw    = S*w_t;
        logSw = log(Sw);
        logSw(Sw == 0) = 0;

        fval(k) = 2*z_k'*w_t - alpha*sum(logSw) + beta*(w_t'*w_t);

        w_opt_k               = w_k;
        w_opt_k(w_opt_k < th) = 0;

        W_opt(:,k) = w_opt_k;

        err(k) = norm(w_opt_k-w_ref);
    
    end

end