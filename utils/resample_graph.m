function Grs = resample_graph(G, model, params, f_samp)

    m = size(G,1);
    num_resample = round(f_samp*nnz(triu(G,1)));

    [rows, cols] = find(triu(G, 1));
    upper_tri_indices = sub2ind([m,m], rows, cols);

    resample_indices = randsample(upper_tri_indices, num_resample);

    for i = 1:num_resample
        idx = resample_indices(i);
        [k,l] = ind2sub([m,m], idx);
        G(k,l) = 0;
        G(l,k) = 0;
    end

    switch lower(model)

        case 'er'
            p = params.var1;

            for i = 1:num_resample
                idx = resample_indices(i);
                [k,l] = ind2sub([m,m], idx);
                G(k,l) = rand < p;
                G(l,k) = G(k,l);
            end

        case 'pa'
            degrees = sum(G, 2);
            prob = degrees/sum(degrees);

            for i = 1:num_resample
                while true
                    k = randsample(1:m, 1, true, prob);
                    l = randsample(1:m, 1, true, prob);
                    if k ~= l && G(k,l) == 0
                        G(k,l) = 1;
                        G(l,k) = 1;
                        break;
                    end
                end
            end

        case 'gaussian'
            s = params.var1;
            T = params.var2;

            corr_mat = normrnd(0, s, m, m);
            corr_mat = (corr_mat + corr_mat')/2;
            corr_mat(1:m+1:end) = 1;

            Grs = abs(corr_mat) > T;
            Grs = triu(Grs,1);
            Grs = Grs + Grs'; 
 
            for i = 1:num_resample
                idx = resample_indices(i);
                [k,l] = ind2sub([m,m], idx);
                G(k,l) = Grs(k,l);
                G(l,k) = G(k,l);
            end

        otherwise
            error('Unknown model type.');
    end

    Grs = G;

end