addpath('./methods/');
addpath('./utils/');
addpath('./configs/');

seed = 123;
rng(seed);

params_table = read_params_table(config_filename);

test_type = params_table.test_type;

if test_type == 0

    graph_model = params_table.graph_model;
    n_samples   = params_table.n_samples;
    n_nodes     = params_table.n_nodes;
    n_blocks    = params_table.n_blocks;
    sigma       = params_table.sigma;

    switch graph_model

        case 'gaussian'
            params_graph.var1 = params_table.var1; % T
            params_graph.var2 = params_table.var2; % S

        case 'er'
            params_graph.var1 = params_table.var1; % p

        case 'pa'
            params_graph.var1 = params_table.var1; % P

    end

    f_resamp = params_table.f_resamp;

else

    data_filename = params_table.data_filename;
    n_samples     = params_table.n_samples;
    sigma         = params_table.sigma;

end

th                = params_table.th;
log10rho_opadmm_i = params_table.log10rho_opadmm_i;
log10rho_opadmm_f = params_table.log10rho_opadmm_f;
n_rho_opadmm      = params_table.n_rho_opadmm;
log10t1_opadmm_i  = params_table.log10t1_opadmm_i;
log10t1_opadmm_f  = params_table.log10t1_opadmm_f;
n_t1_opadmm       = params_table.n_t1_opadmm;
log10t2_opadmm_i  = params_table.log10t2_opadmm_i;
log10t2_opadmm_f  = params_table.log10t2_opadmm_f;
n_t2_opadmm       = params_table.n_t2_opadmm;
pwr_vec_i         = params_table.pwr_vec_i;
pwr_vec_f         = params_table.pwr_vec_f;
pwr_vec_step      = params_table.pwr_vec_step;
y_range_i         = params_table.y_range_i;
y_range_f         = params_table.y_range_f;

pwr_vec = pwr_vec_i:pwr_vec_step:pwr_vec_f;

params_ref.tol       = 1e-12;
params_ref.maxit     = 10000;
params_ref.step_size = 0.2;

y_range = [y_range_i,y_range_f];
fig_idx = 1;

params_signal.n_samples = n_samples;
params_signal.sigma     = sigma;

params_gsearch.n_samples = n_samples;
params_gsearch.th        = th;
params_gsearch.pwr_vec   = pwr_vec;

rho_opadmm_v = 10.^(log10rho_opadmm_i:(log10rho_opadmm_f-log10rho_opadmm_i)/(n_rho_opadmm-1):log10rho_opadmm_f);
t1_opadmm_v  = 10.^(log10t1_opadmm_i:(log10t1_opadmm_f-log10t1_opadmm_i)/(n_t1_opadmm-1):log10t1_opadmm_f);
t2_opadmm_v  = 10.^(log10t2_opadmm_i:(log10t2_opadmm_f-log10t2_opadmm_i)/(n_t2_opadmm-1):log10t2_opadmm_f);

err_opadmm_m = zeros([n_samples,n_rho_opadmm,n_t2_opadmm,n_t1_opadmm]);

color_v = generate_contrast_colors(max([n_rho_opadmm,n_t2_opadmm,n_t1_opadmm]));
ticks   = n_samples/5:n_samples/5:n_samples;

params_signal.n_samples = n_samples;
params_signal.sigma     = sigma;

%%

if test_type == 0

    params_graph.n_nodes = n_nodes;
    [G,X] = gen_graphsignals(graph_model, f_resamp, n_blocks, params_signal, params_graph);
    [alpha,beta,W_opt] = graphlearn_gridsearch(G, X, n_blocks, params_gsearch, params_ref);

else

    [G,X] = gen_graphsignals_from_data(params_signal, data_filename);
    [alpha,beta,W_opt] = graphlearn_gridsearch(G, X, 1, params_gsearch, params_ref);
    n_nodes = size(G{1},1);

end

[S,St] = sum_squareform(n_nodes);
maxeval_StS = 2*(n_nodes-1);

%%

clear p;

param_solver.gamma       = 2e-3;
param_solver.inner_loops = 1;
param_solver.n_samp      = 1;
param_solver.th          = 1e-5;

for k = 1:n_rho_opadmm

    param_solver.rho_0 = rho_opadmm_v(k);

    tau1_opadmm_str = cellfun(@(x) sprintf('t_1=%g',x), num2cell(t1_opadmm_v), 'UniformOutput', false);
    tau2_opadmm_str = cellfun(@(x) sprintf('t_2=%g',x), num2cell(t2_opadmm_v), 'UniformOutput', false);

    for j = 1:n_t2_opadmm

        param_solver.tau2_fun = @(rho) t2_opadmm_v(j)/rho;

        for l = 1:n_t1_opadmm

            param_solver.tau1_fun = @(rho) t1_opadmm_v(l)/(rho*maxeval_StS);
            [~,~,err] = opadmm_gl(X, alpha, beta, W_opt, param_solver);

            figure(fig_idx); subplot(1,n_t2_opadmm,j); hold on;
            p(l) = plot(err,'-','LineWidth',1.0,'Color',color_v(l,:)); hold off;

            err_opadmm_m(:,k,j,l) = err;

        end

        grid on;
        ax = gca;
        ax.GridLineStyle = '--';
        ax.GridAlpha = 0.75;
        xticks(ticks);
        ylim(y_range);
        legend(p, tau1_opadmm_str, 'Location', 'southwest');
        xlabel('Number of iterations');
        ylabel('$\|\mathbf{w} - \hat{\mathbf{w}}\|_{2}$', 'Interpreter', 'Latex');
        set(gca, 'YScale', 'log');
        title(sprintf('t_2=%g', t2_opadmm_v(j)));

    end

    sgtitle(sprintf('\\rho=%g', rho_opadmm_v(k)));

    fig_idx = fig_idx+1;

end

%% comparison and grid search

clear p;

% parameters for all solvers
param_solver.gamma         = 2e-3;
param_solver.inner_loops   = 1;
param_solver.n_samp        = 1;
param_solver.th            = 1e-5;

% online PG
[~,~,err_opg] = opg_gl(X, alpha, beta, W_opt, param_solver);

% online DPG
[~,~,err_odpg] = odpg_gl(X, alpha, beta, W_opt, param_solver);

% online P-ADMM (grid search)
err_ref_m = repmat(err_odpg, 1, n_rho_opadmm, n_t2_opadmm, n_t1_opadmm);
diff_err  = err_ref_m - err_opadmm_m;
diff_err  = diff_err(round(0.9*n_samples):end,:,:,:);
diff_posi = squeeze(diff_err(end,:,:,:) > 0);
gs_val    = squeeze(sum(diff_err.*(diff_err > 0), 1));

if sum(diff_posi(:)) > 0
    gs_val = diff_posi.*gs_val;
end

[~,idx_max] = max(gs_val(:));
[k_max,j_max,l_max] = ind2sub(size(gs_val), idx_max);

err_opadmm = err_opadmm_m(:,k_max,j_max,l_max);

rho_opadmm_opt = rho_opadmm_v(k_max);
t2_opadmm_opt  = t2_opadmm_v(j_max);
t1_opadmm_opt  = t1_opadmm_v(l_max);

%% plots

h1 = figure(fig_idx);
hold on;
p1 = plot(err_opg, '-', 'LineWidth', 2.5,'Color', 'r');
hold on;
p2 = plot(err_odpg, '-', 'LineWidth', 2.5,'Color', 'b');
hold on;
p3 = plot(err_opadmm, '-', 'LineWidth', 2.5,'Color', 'm');

h1.Position = [100, 100, 560, 455];

axes_FontSize = 18;
ticks_FontSize = 16;

grid on;
ax = gca;
box(ax,'on');

ax.GridLineStyle = '--';
ax.GridAlpha = 0.75;

legend([p1 p2 p3], {'online PG','online DPG','OPADMM'},'Location','southwest','FontSize',axes_FontSize);
xlabel('$k$','Interpreter','Latex','FontSize',axes_FontSize);
ylabel('$\|\mathbf{w}^{(k)} - \hat{\mathbf{w}}\|_{2}$','Interpreter','Latex','FontSize',axes_FontSize);

xticks(ticks);

xlim([1,n_samples]);
ylim(y_range);

set(gca, 'YScale','log');

x_axis = get(gca,'XTickLabel');
y_axis = get(gca,'YTickLabel');
set(gca,'XTickLabel',x_axis,'fontsize',ticks_FontSize);
set(gca,'YTickLabel',y_axis,'fontsize',ticks_FontSize);