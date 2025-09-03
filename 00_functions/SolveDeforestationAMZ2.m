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

    Candidates = Dnext & ~Abiz & border ==1;  % plots that can gain rights

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

estados = readgeotable("..\02_inputs\states_legal_amazon.shp");

%% ---- INPUTS expected ----
% D(:,:,1)     : initial deforestation (0/1)
% Dnow         : final deforestation (0/1)
% Protected    : protected raster (0/1), used for both unless Pnow provided
% (optional) Pnow : protected raster at final time (0/1)
% Agribusiness : initial agribusiness raster (0/1)
% Abiz         : final agribusiness raster (0/1)
% mask         : raster mask (1 = on map, 0 = off map)
% estados      : vector layer with state borders (e.g., shaperead mapstruct)

%% ---- Derive masks & sanity checks ----
%% Expected inputs:
% D(:,:,1)  -> initial deforestation (0/1)
% Dnow      -> final deforestation (0/1)
% Protected -> protected areas (0/1)  [optionally Pnow; otherwise Protected reused]
% Agribusiness -> initial agribusiness (0/1)
% Abiz         -> final agribusiness (0/1)
% mask      -> 1 inside map, 0 outside (same size as D(:,:,1))
% estados   -> vector borders (mapstruct with X/Y or Lon/Lat, polyshape, or geotable-like)

%% Masks (logical) and consistency
D0 = logical(D(1:params.m,1:params.n,1));   D1 = logical(Dnow(1:params.m,1:params.n));
P0 = logical(L(1:params.m,1:params.n));
P1 = P0;
A0 = Agribusiness(1:params.m,1:params.n,1);
A1 = A0;

F0 = ~D0;   F1 = ~D1;

mask01 = logical(mask(1:params.m,1:params.n));

% Erase off-map
F0 = F0 & mask01;  P0 = P0 & mask01;  D0 = D0 & mask01;  A0 = A0 & mask01;
F1 = F1 & mask01;  P1 = P1 & mask01;  D1 = D1 & mask01;  A1 = A1 & mask01;

[m,n] = size(D0);

%% Colors (Forest, Deforested, Protected, Agribusiness)
clrForest = [0.00 0.55 0.00];   % darker green
clrDefor  = [0.55 0.20 0.75];   % (unchanged) red
clrProt   = [0.50 0.72 1.00];   % lighter blue
clrAgri   = [1.00 0.00 0.00];   % purple (replaces yellow)
greyLines = [0.15 0.15 0.15];   % darker neutral grey for estados lines (optional)
alphaLay  = 0.95;               % (unchanged)

% Painter for a solid truecolor layer with alpha mask
paintLayer = @(ax,mask3,color,alphaVal) image(ax, ...
    'CData', repmat(reshape(color,1,1,3), [m n 1]), ...
    'AlphaData', double(mask3) * alphaVal);

%% Figure with two maps and shared legend at bottom
t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% ---------- Left: Initial ----------
ax1 = nexttile(t,1);
hold(ax1,'on'); axis(ax1,'image'); set(ax1,'YDir','reverse');  % correct orientation
title(ax1,'Initial Deforestation');
% base white (so masked-out is white)
image(ax1,'CData',ones(m,n,3),'AlphaData',ones(m,n));

% Draw in required order: Forest -> Protected -> Deforested -> Agribusiness
paintLayer(ax1, F0, clrForest, alphaLay);
paintLayer(ax1, P0, clrProt,   alphaLay);
paintLayer(ax1, D0, clrDefor,  alphaLay);
paintLayer(ax1, A0, clrAgri,   alphaLay);

% estados on top (darker grey)
plot_estados_top(ax1, estados, greyLines);

axis(ax1,'tight'); ax1.XTick=[]; ax1.YTick=[]; box(ax1,'on'); hold(ax1,'off');

% ---------- Right: Final ----------
ax2 = nexttile(t,2);
hold(ax2,'on'); axis(ax2,'image'); set(ax2,'YDir','reverse');
title(ax2,'Final Deforestation');
image(ax2,'CData',ones(m,n,3),'AlphaData',ones(m,n));

paintLayer(ax2, F1, clrForest, alphaLay);
paintLayer(ax2, P1, clrProt,   alphaLay);
paintLayer(ax2, D1, clrDefor,  alphaLay);
paintLayer(ax2, A1, clrAgri,   alphaLay);

plot_estados_top(ax2, estados, greyLines);

axis(ax2,'tight'); ax2.XTick=[]; ax2.YTick=[]; box(ax2,'on'); hold(ax2,'off');

% ---------- Shared bottom legend (force 4 entries) ----------
hold(ax2,'on');  % create legend handles here
hL(1) = plot(ax2, NaN,NaN,'s','MarkerFaceColor',clrForest,'MarkerEdgeColor','k','LineStyle','none','MarkerSize',10,'DisplayName','Forest');
hL(2) = plot(ax2, NaN,NaN,'s','MarkerFaceColor',clrProt  ,'MarkerEdgeColor','k','LineStyle','none','MarkerSize',10,'DisplayName','Protected');
hL(3) = plot(ax2, NaN,NaN,'s','MarkerFaceColor',clrDefor ,'MarkerEdgeColor','k','LineStyle','none','MarkerSize',10,'DisplayName','Deforested');
hL(4) = plot(ax2, NaN,NaN,'s','MarkerFaceColor',clrAgri  ,'MarkerEdgeColor','k','LineStyle','none','MarkerSize',10,'DisplayName','Agribusiness');
hold(ax2,'off');

lg = legend(hL, 'Orientation','horizontal','NumColumns',4,'AutoUpdate','off');
lg.Layout.Tile = 'south';


%% -------- Deforestation & Emissions Time Series --------
figure;
subplot(1,2,1)
plot(1:params.Tmax-1, deforest(1:params.Tmax-1), '-o', 'LineWidth', 1.5);
xlabel('Time'); ylabel('Cumulative Deforestation');
title('Deforestation over Time'); grid on;

subplot(1,2,2)
plot(1:params.Tmax-1, emissions(1:params.Tmax-1), '-o', 'LineWidth', 1.5);
xlabel('Time'); ylabel('Cumulative Emissions');
title('Cumulative Emissions over Time'); grid on;

ret = D;
