function [W_opt, fval, err] = odpg_gl(X, alpha, beta, W_ref, param_solver)

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

    z_k   = zeros(p,1);
    lmb_k = zeros(m,1);

    W_opt = zeros(p,n);

    fval = zeros(n,1);
    err  = zeros(n,1);

    [S,St] = sum_squareform(m);

    L  = (m-1)/beta;
    mu = 1/L;

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

            v_k = max(0,(St*lmb_k-2*z_k)/(2*beta));
 
            d_k = S*v_k - L*lmb_k;
            u_k = (d_k + sqrt(d_k.^2+4*alpha*L))/2;
            
            lmb_km1 = lmb_k;
            lmb_k   = lmb_km1 - mu*(S*v_k-u_k);

        end

        fprintf('Iteration %d completed! Elapsed time: %.2f sec.\n', k, toc(t_0));

        w_k = max(0,(St*lmb_k-2*z_k)/(2*beta));

        Sw             = S*w_k;
        logSw          = log(Sw);
        logSw(Sw == 0) = 0;

        fval(k) = 2*z_k'*w_k - alpha*sum(logSw) + beta*(w_k'*w_k);

        w_opt_k               = w_k;
        w_opt_k(w_opt_k < th) = 0;

        W_opt(:,k) = w_opt_k;

        err(k) = norm(w_opt_k-w_ref);

    end

end