clear;

%% 1) 正向模型与真值
sigma_ground = 0.6;
sigma_background = 0.12;
drive_electrodes = [601 609 617 625 425 225 25 17 9 1 201 401];

fwd_mdl = mk_common_model('d2s', 12);
fwd_mdl.fwd_model = assign_electrode_nodes(fwd_mdl.fwd_model, drive_electrodes);

img_bg = mk_image(fwd_mdl.fwd_model, sigma_background);
img_true = apply_conductivity_pattern(img_bg, sigma_ground, sigma_background);

vh = fwd_solve(img_bg);
vi = fwd_solve(img_true);

%% 2) 反演模板
measurement_electrodes = [9313 9345 9377 9409 6208 3201 97 65 33 1 3105 6209];
inv_base = mk_common_model('i2s', 12);
inv_base.fwd_model = assign_electrode_nodes(inv_base.fwd_model, measurement_electrodes);

inv_template = eidors_obj('inv_model', 'Normalized absolute reconstruction');
inv_template.reconst_type = 'absolute';
inv_template.jacobian_bkgnd.value = 1;
inv_template.fwd_model = inv_base.fwd_model;
inv_template.solve = @inv_solve_core;
inv_template.inv_solve_core.max_iterations = 10;
inv_template.inv_solve_core.term_tolerance = 1e-9;
inv_template.inv_solve_core.verbose = 1;

%% 3) 先做标准TV作为 pilot（统一超参数，保证公平对比）
hp = 0.05; % 与下文的形态路径先验保持一致
imdl_tv = inv_template;
imdl_tv.hyperparameter = struct('value', hp);
imdl_tv.RtR_prior = @(m) prior_TV(m)' * prior_TV(m);
if isfield(imdl_tv, 'R_prior'); imdl_tv = rmfield(imdl_tv, 'R_prior'); end
imgr_tv = inv_solve(imdl_tv, vi);

%% 4) 三组形态路径先验参数（[thr_factor, reduce_factor, morph_thresh]）
param_sets = [
    0.3, 0.5, 0.2;  % 更宽松的阈值 + 路径较弱衰减
    0.5, 0.2, 0.3;  % 基准（与 run_morph_path_tv 中一致）
    0.7, 0.2, 0.4   % 更严格阈值，路径更稀
];

results = cell(size(param_sets,1),1);
labels  = strings(size(param_sets,1),1);

for k = 1:size(param_sets,1)
    thr = param_sets(k,1); red = param_sets(k,2); mth = param_sets(k,3);

    imdl_morph = inv_template;
    imdl_morph.hyperparameter = struct('value', hp); % 与 pilot 相同的正则超参数
    % 用 pilot 结果引导路径检测
    imdl_morph.elem_data = imgr_tv.elem_data;
    imdl_morph.RtR_prior = @(inv_model) build_path_regularization_RtR_with_morph(inv_model, thr, red, mth);
    if isfield(imdl_morph, 'R_prior'); imdl_morph = rmfield(imdl_morph, 'R_prior'); end

    % --- Diagnostics: call builder directly to inspect mask/weights/RtR ---
    try
        [RtRtest, diagnostics] = build_path_regularization_RtR_with_morph(imdl_morph, thr, red, mth);
        fprintf('Param set %d: thr=%.2f reduce=%.2f morph=%.2f -> mask2=%d, path_nodes=%d\n', ...
            k, thr, red, mth, sum(diagnostics.mask2), numel(diagnostics.path_nodes));
        % robust condition estimate for sparse RtR
        try
            if issparse(RtRtest)
                kappa_est = condest(RtRtest); % 1-norm condition number estimate
                rcond_est = 1./kappa_est;
            else
                rcond_est = rcond(RtRtest);
            end
        catch
            rcond_est = NaN;
        end
        fprintf('  w: min/med/max = %.4g / %.4g / %.4g, RtR rcond_est~=%.3e\n', min(diagnostics.w), median(diagnostics.w), max(diagnostics.w), rcond_est);

        % show small diagnostic figures (non-blocking)
        figure('Name',sprintf('diag_g_elem2_set%d',k)); show_fem(mk_image(imdl_morph.fwd_model, diagnostics.g_elem2)); title(sprintf('g_elem2 (set %d)',k));
        figure('Name',sprintf('diag_mask2_set%d',k)); show_fem(mk_image(imdl_morph.fwd_model, double(diagnostics.mask2))); title(sprintf('mask2 (set %d)',k));
        % map edge weights to element-averaged weights for display
        EInc = spones(abs(prior_TV(imdl_morph)));
        edges_per_elem = full(max(1, sum(EInc,1)'));
        elem_weight = (EInc' * diagnostics.w) ./ edges_per_elem;
        figure('Name',sprintf('diag_elemw_set%d',k)); show_fem(mk_image(imdl_morph.fwd_model, elem_weight)); title(sprintf('elem weight (set %d)',k));
    catch ME
        warning('Diagnostics failed for set %d: %s', k, ME.message);
        RtRtest = [];
        diagnostics = [];
    end

    % Ensure cache won't return a previous run: clear eidors cache if present
    try
        if exist('eidors_cache','file'), eidors_cache('clear_all'); end
    catch
    end

    % Now run the actual inversion for this parameter set
    results{k} = inv_solve(imdl_morph, vi);
    labels(k) = sprintf('thr=%.2f, reduce=%.2f, morph=%.2f', thr, red, mth);
end

%% 5) 可视化
figure('Name','Morphological Path-TV: parameter sweep');
subplot(2,2,1); show_slices(imgr_tv); title('Standard TV (pilot)');
for k = 1:size(param_sets,1)
    subplot(2,2,k+1); show_slices(results{k}); title(labels(k));
end

%% 6) 可选：打印每个结果的测量残差
try
    pred_tv = fwd_solve(imgr_tv); 
    tv_vec = get_meas_vector(pred_tv);
    vi_vec = get_meas_vector(vi);
    fprintf('Std TV residual (L2): %.6g (rel=%.6g)\n', norm(tv_vec-vi_vec), norm(tv_vec-vi_vec)/norm(vi_vec));
    for k = 1:numel(results)
        pk = fwd_solve(results{k}); pv = get_meas_vector(pk);
        fprintf('Morph %d residual (L2): %.6g (rel=%.6g)\n', k, norm(pv-vi_vec), norm(pv-vi_vec)/norm(vi_vec));
    end
catch
end

function v = get_meas_vector(s)
if isstruct(s)
    if isfield(s,'meas'), v = s.meas; elseif isfield(s,'v'), v = s.v; elseif isfield(s,'volt'), v = s.volt; else
        c = struct2cell(s); idx = find(cellfun(@isnumeric,c)); v = c{idx(1)}; 
    end
else
    v = s;
end
end
