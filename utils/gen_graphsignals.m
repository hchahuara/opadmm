function [G,X] = gen_graphsignals(graph_model, f_resamp, n_blocks, params_signal, params_graph)

    n_samples = params_signal.n_samples;
    sigma     = params_signal.sigma;

    n_nodes   = params_graph.n_nodes;

    switch graph_model

        case 'gaussian'
            var1 = params_graph.var1;
            var2 = params_graph.var2;
            [G_k,~,~] = construct_graph(n_nodes, graph_model, var1, var2);

        case 'er'
            var1 = params_graph.var1;
            [G_k,~,~] = construct_graph(n_nodes, graph_model, var1);

        case 'pa'
            var1 = params_graph.var1;
            [G_k,~,~] = construct_graph(n_nodes, graph_model, var1);

        case 'ff'
            var1 = params_graph.var1;
            var2 = params_graph.var2;
            [G_k,~,~] = construct_graph(n_nodes, graph_model, var1, var2);

        case 'chain'
            [G_k,~,~] = construct_graph(n_nodes, graph_model);

    end

    L = diag(sum(full(G_k),2)) - G_k; % laplacian

    [V,D] = eig(L);

    X_gft = mvnrnd(zeros(1,n_nodes), pinv(D), n_samples/n_blocks); % random gft coefficients
    X_teo = V*X_gft'; % signal on graph G
    X_k   = X_teo + sigma*randn(size(X_teo)); % noisy signal

    G{1} = G_k;
    X    = X_k;

    for k = 2:n_blocks

        G_km1 = G_k;
        G_k    = resample_graph(G_km1, graph_model, params_graph, f_resamp);

        L = diag(sum(full(G_k),2)) - G_k; % laplacian

        [V,D] = eig(L);

        X_gft = mvnrnd(zeros(1,n_nodes), pinv(D), n_samples/n_blocks); % random gft coefficients
        X_teo = V*X_gft'; % signal on graph G
        X_k   = X_teo + sigma*randn(size(X_teo)); % noisy signal

        G{k} = G_k;
        X    = cat(2, X, X_k);

    end

end