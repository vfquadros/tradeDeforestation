function ret = SolveDeforestationAMZ(params,A,F,O,D,L,Agribusiness,O_aux,w1,animation,mask)

q1 = (1-params.beta.*(1-params.pi))./(params.pi.*(params.beta.^params.Tlag));
q2 = (1-params.beta.*(1-params.pi).*(1+params.pi.*(params.beta.^params.Tlag)))./(params.pi.*(params.beta.^params.Tlag).*(1-params.beta));
q3 = (params.alpha.*(1-params.beta.^params.Tlag)+(params.beta.^params.Tlag).*(1-params.pi).*(1-params.beta))./(params.pi.*(params.beta.^params.Tlag).*(1-params.beta));

deforest = zeros(1,params.Tmax);
emissions = zeros(1,params.Tmax);
DiffProfit = zeros(1,params.N);  % This will store the differences in Agribusinessess

if animation == true% function psi in the model
    figure;                % Starting animation
    h = imagesc(D(:,:,1)); % Create heatmap
    colorbar;
    set(h, 'CData', D(:,:,1)); % Update heatmap data
    pause(1.5); % Pause for smooth animation
end
time = tic();
for t = 1:params.Tmax-1
    % ===== Step 1: calculating agribusiness diff profits =====
    B = (Agribusiness(:,:,t) == 1);                % current agribusiness map (logical)
    
    % -- find 8-connected clusters (requires Image Processing Toolbox) --
    CC = bwconncomp(B, 8);
    clusters = CC.PixelIdxList;                    % cell array of linear indices
    
    % -- border neighbors of each cluster via convolution (4-neighborhood) --
    K8 = ones(3); K8(2,2)=0;                    % Von Neumann neighborhood
    neighbors = cell(size(clusters));
    for i = 1:numel(clusters)
        M = false(params.N,params.N); M(clusters{i}) = true;
        borderMask = conv2(double(M), K8, 'same') > 0 & ~M;   % adjacent but not in cluster
        neighbors{i} = find(borderMask);                      % linear indices
    end
    
    % -- compute DiffProfit only on border cells not already agribusiness --
    DiffProfit = zeros(params.N,params.N);
    for i = 1:numel(neighbors)
        clusterIdx = clusters{i};
        ContAgri   = numel(clusterIdx);                       % size of cluster
    
        % profit() should accept linear indices; term1/term2 are scalars
        term1 = sum(profit(params,A, ContAgri+1, clusterIdx));
        term2 = sum(profit(params,A, ContAgri,   clusterIdx));
    
        nb = neighbors{i};
        nb = nb(~B(nb));                                      % exclude cells already in agribusiness
    
        % fill marginal diff profit on border candidates
        DiffProfit(nb) = term1 - term2 + profit(params,A, ContAgri+1, nb);
    end
    
    % ===== Step 2 (NxN): farmers' profit and price R =====
    FarmersProfit = max(q1.*F + q2.*params.wage - q3.*A - params.wage/(1-params.beta), (A - params.wage)/(1-params.beta));
    R(:,:,t) = w1.*DiffProfit + (1-w1).*FarmersProfit;              % adjust weights as desired
    %q2
    %DiffProfit(27:33,24:33)
    %FarmersProfit(117:123,117:123)
    % ===== Farmers' decisions (NxN) =====
    Vu = value_u(params,A, R(:,:,t), F);
    Vp = value_p(params,A, F);
    
    Du = (L == 0) & (D(:,:,t) ~= 1) & (Agribusiness(:,:,t) ~= 1) & (Vu > 0);
    Dp = (L == 1) & (D(:,:,t) ~= 1) & (Agribusiness(:,:,t) ~= 1) & (Vp > 0);
    
    % ===== State updates (NxN) =====
    D(:,:,t+1) = D(:,:,t) + Du + Dp;
    O(:,:,t+1) = O(:,:,t) + D(:,:,t) .* O_aux;               % O_aux must be N x N
    
    Candidates = D(:,:,t) - Agribusiness(:,:,t);              % plots that can gain rights
    PP = (Candidates == 1) & (rand(params.N,params.N) < params.pi) & (L == 0);      % earn property rights
    NewAgri = (PP == 1) & (DiffProfit >= R(:,:,t)) & (mask==1);
    Agribusiness(:,:,t+1) = Agribusiness(:,:,t) + NewAgri;
    
    % (optional) visualization + aggregates
    if animation == true
        set(h, 'CData', D(:,:,t));                                % heatmap update
        pause(0.1); % Pause for smooth animation
    end
    deforest(t)  = sum(D(:,:,t), 'all');
    emissions(t) = sum(O(:,:,t), 'all');
end
toc(time)

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
imagesc(D(:,:,end)); colormap([0 1 0; 1 0 0]);   % green for 0, red for 1
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
alphaMask = double(L == 1 & squeeze(D(:,:,end)) == 0);

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

ret = true;