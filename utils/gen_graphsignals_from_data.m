function [G,X] = gen_graphsignals_from_data(params_signal,data_filename)

    n_samples = params_signal.n_samples;
    sigma     = params_signal.sigma;

    data = load(data_filename);
    %data = load('./data/bcspwr01.mat');
    %data = load('./data/mesh1e1.mat');
    %data = load('./data/bcspwr03.mat');
    %data = load('./data/lshp_265.mat');

    G_k = data.Problem.A;

    n_nodes = size(G_k,1);

    L = diag(sum(full(G_k),2)) - G_k; % laplacian

    [V,D] = eig(L);

    X_gft = mvnrnd(zeros(1,n_nodes), pinv(D), n_samples); % random gft coefficients
    X_teo = V*X_gft'; % signal on graph G
    X_k   = X_teo + sigma*randn(size(X_teo)); % noisy signal

    G{1} = G_k;
    X    = X_k;

end