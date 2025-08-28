function ret = SolveDeforestationAMZ2(params,A,F,O,D,L,Agribusiness,O_aux,w1,animation,mask)

% Precompute some functions
q1 = (1-params.beta.*(1-params.pi))./(params.pi.*(params.beta.^params.Tlag));
q2 = (1-params.beta.*(1-params.pi).*(1+params.pi.*(params.beta.^params.Tlag)))./(params.pi.*(params.beta.^params.Tlag).*(1-params.beta));
q3 = (params.alpha.*(1-params.beta.^params.Tlag)+(params.beta.^params.Tlag).*(1-params.pi).*(1-params.beta))./(params.pi.*(params.beta.^params.Tlag).*(1-params.beta));

% Assume: D, O, Agribusiness, A, F, L, mask are N×N (logical or single where possible)
% Keep rolling buffers to save RAM:
Dnow  = D(:,:,1);    Onow = O(:,:,1);
Abiz  = Agribusiness(:,:,1);
% Optional: keep Rnow only (no full R history)
% Precompute constants
c_wage = params.wage/(1-params.beta);

deforest  = zeros(params.Tmax-1,1);
emissions = zeros(params.Tmax-1,1);

se = strel('square',3);

if animation == true% function psi in the model
    figure;                % Starting animation
    h = imagesc(D(:,:,1)); % Create heatmap
    colorbar;
    set(h, 'CData', D(:,:,1)); % Update heatmap data
    pause(1.5); % Pause for smooth animation
end
time = tic();
for t = 1:params.Tmax-1
    % ===== Step 1: agribusiness marginal diff profits =====
    B = Abiz; % current agribusiness map (logical)

    % Connected components (8-neighborhood) and labels
    CC  = bwconncomp(B,8);
    LCC = labelmatrix(CC);    % int32 labels, 0 where B==0

    % Border of agribusiness (adjacent but not in B)
    border = imdilate(B,se) & ~B;   % logical N×N

    DiffProfit = zeros(size(B),'like',A);  % single/double; match A type

    % Prepare parallel accumulation (option 1: cell of sparse)
    dpCells = cell(CC.NumObjects,1);

    parfor k = 1:CC.NumObjects
        pix  = CC.PixelIdxList{k};        % linear indices of cluster k
        Sk   = numel(pix);                % cluster size

        % Compute cluster terms once
        % Vectorized profit evals; profit returns vector for indices
        p1 = profit(params, A, Sk+1, pix);
        p0 = profit(params, A, Sk,   pix);
        term1_minus_term2 = sum(p1) - sum(p0);  % scalar

        % Border neighbors mask for THIS cluster (no neighbor list!)
        Mk = false(size(B)); Mk(pix) = true;
        nbMask = (imdilate(Mk, se) & ~Mk) & border;  % pixels bordering this cluster

        if any(nbMask(:))
            nb = find(nbMask);

            % Profit term on border candidates with S=Sk+1 (vectorized)
            p_nb = profit(params, A, Sk+1, nb);  % vector

            % Local buffer (sparse is memory-lean for scattered writes)
            local = sparse(nb, ones(numel(nb),1), term1_minus_term2 + p_nb, numel(B), 1);

            dpCells{k} = local;
        else
            dpCells{k} = sparse(numel(B),1);
        end
    end

    % Combine partials once
    if ~isempty(dpCells)
        DiffProfit(:) = full(sum(cat(2, dpCells{:}), 2));
    % (optional) visualization + aggregates
    if animation == true
        set(h, 'CData', Dnow);                                % heatmap update
        pause(0.1); % Pause for smooth animation
    end
    end
    toc(time)

    % ===== Step 2: Farmers' profit & R (no 3D R kept) =====
    % As vectors/matrices; ensure elementwise ops
    FarmersProfit = max( q1.*F + q2.*params.wage - q3.*A - c_wage, ...
                         (A - params.wage)./(1-params.beta) );

    Rnow = w1.*DiffProfit + (1-w1).*FarmersProfit;

    % ===== Farmers' decisions =====
    Vu = value_u(params, A, Rnow, F);   % vectorized inside
    Vp = value_p(params, A, F);

    Du = (L==0) & ~Dnow & (Vu > 0) & mask==1;
    Dp = (L==1) & ~Dnow & (Vp > 0) & mask==1;

    % ===== State updates =====
    Dnext = Dnow | Du | Dp;

    % O update (use in-place math; O_aux is N×N)
    aux = O_aux;
    aux(~Dnow) = 0;
    Onext = Onow + aux;

    Candidates = Dnext & ~Abiz;  % plots that can gain rights

    % Random only on candidate set (saves time & memory)
    candIdx = find(Candidates & (L==0));
    keep    = false(size(Dnow));
    if ~isempty(candIdx)
        keep(candIdx(rand(numel(candIdx),1) < params.pi)) = true; % PP
    end

    PP = keep;   % logical
    NewAgri = PP & (DiffProfit >= Rnow) & (mask==1);

    Abiz_next = Abiz | NewAgri;

    % (optional) animation
    if animation
        set(h,'CData', Dnow); drawnow limitrate;
    end

    % Aggregates (ensure no NaNs if possible)
    deforest(t)  = sum(Dnow,'all');
    emissions(t) = sum(Onext,'all');

    % Roll
    Dnow  = Dnext;
    Onow  = Onext;
    Abiz  = Abiz_next;

    % If you must keep history:
    % D(:,:,t+1) = Dnext; O(:,:,t+1) = Onext; Agribusiness(:,:,t+1) = Abiz_next;
end

% === Final Heatmaps ===
% Figure 1: Initial x Final Feforestiation
figure;
subplot(1,2,1)

% --- Base map (D) ---
imagesc(D(:,:,1));
colormap([0 1 0; 1 0 0]);   % green for 0, red for 1
clim([0 1]);
hold on;

% Remove continuous colorbar
colorbar('off');

% Figure 2
% --- Overlay (L) only where D==0 ---
overlayRGB = zeros([size(L), 3]);  % initialize RGB image
overlayRGB(:,:,3) = 1;             % blue channel = 1

% Alpha mask: 1 only where L==1 and D==0
alphaMask = double(L == 1 & squeeze(D(:,:,1)) == 0);

% Expand alpha mask to match image size if needed
hOverlay = imagesc(overlayRGB);
set(hOverlay, 'AlphaData', alphaMask);

hold off;
title('Initial Deforestation');

subplot(1,2,2)
imagesc(Dnow); colormap([0 1 0; 1 0 0]);   % green for 0, red for 1
clim([0 1]);
hold on;

% Remove continuous colorbar
colorbar('off');

% Add a custom discrete legend
h = zeros(3,1); % placeholder for legend handles
h(1) = plot(NaN,NaN,'s','MarkerFaceColor',[0 1 0],'MarkerEdgeColor','k'); % green square
h(2) = plot(NaN,NaN,'s','MarkerFaceColor',[1 0 0],'MarkerEdgeColor','k'); % red square
h(3) = plot(NaN,NaN,'s','MarkerFaceColor',[0 0 1],'MarkerEdgeColor','k'); % blue square
legend(h,{'Forest','Deforested', 'Protected'},'Location','bestoutside');

% --- Overlay (L) only where D==0 ---
overlayRGB = zeros([size(L), 3]);  % initialize RGB image
overlayRGB(:,:,3) = 1;             % blue channel = 1

% Alpha mask: 1 only where L==1 and D==0
alphaMask = double(L == 1 & squeeze(Dnow) == 0);

% Expand alpha mask to match image size if needed
hOverlay = imagesc(overlayRGB);
set(hOverlay, 'AlphaData', alphaMask);

hold off;
title('Final Deforestation');

figure;
subplot(1,2,1)
plot(1:params.Tmax-1, deforest(1:params.Tmax-1), '-o', 'LineWidth', 1.5);
xlabel('Time');
ylabel('Cumulative Deforestation');
title('Deforestation over Time');
grid on;

subplot(1,2,2)
plot(1:params.Tmax-1, emissions(1:params.Tmax-1), '-o', 'LineWidth', 1.5);
xlabel('Time');
ylabel('Cumulative Emissions');
title('Cumulative Emissions over Time');
grid on;

ret = D;