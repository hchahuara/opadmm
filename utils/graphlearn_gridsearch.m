function [alpha_opt,beta_opt,W_opt] = graphlearn_gridsearch(G, X, n_blocks, params_gsearch, params_ref)

    n_samples = params_gsearch.n_samples;
    th        = params_gsearch.th;
    pwr_vec   = params_gsearch.pwr_vec;

    F         = zeros(length(pwr_vec)^2, 1);
    hyp_param = zeros(2, length(pwr_vec)^2);

    fprintf('Begining grid search... \n');

    n_samp_per_block = n_samples/n_blocks;

    Z = sparse(gsp_distanz(X(:,1:n_samp_per_block)').^2)/n_samp_per_block;

    alpha = 10.^pwr_vec;
    beta  = 10.^pwr_vec;

    idx = 1;

    for k = 1:length(alpha)
        for l = 1:length(beta)
            hyp_param(:,idx) = [alpha(k);beta(l)];
            idx = idx+1;
        end
    end

    for k = 1:size(hyp_param,2)
        [W_test,~] = gsp_learn_graph_log_degrees(Z, hyp_param(1,k), hyp_param(2,k));
        W_test(W_test < th) = 0;
        [~,~,F(k),~,~] = graph_learning_perf_eval(G{1}, W_test);
    end

    idx    = find(F == max(F));
    temp_b = hyp_param(2,idx);
    b_id   = find(temp_b == max(temp_b));
    idx    = idx(b_id(1));

    alpha_opt = hyp_param(1,idx);
    beta_opt  = hyp_param(2,idx);

    fprintf('The grid search is finished! \n');

    for k = 1:n_blocks
        Z = sparse(gsp_distanz(X(:,(k-1)*n_samp_per_block+1:k*n_samp_per_block)').^2)/n_samp_per_block;
        [W,~] = gsp_learn_graph_log_degrees(Z, alpha_opt, beta_opt, params_ref);
        W(W < th) = 0;
        Wt = squareform_sp(W)';
        W_opt{k} = Wt(:);
    end

end