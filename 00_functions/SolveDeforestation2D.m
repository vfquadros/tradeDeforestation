function SolveDeforestation2D(N,Tmax,A,F,O,D,L,Agribusiness,O_aux,w1)

global N beta alpha pi wage Tlag

q1 = (1-beta.*(1-pi))./(pi.*(beta.^Tlag));
q2 = (1-beta.*(1-pi).*(1+pi.*(beta.^Tlag)))./(pi.*(beta.^Tlag).*(1-beta));
q3 = (alpha.*(1-beta.^Tlag)+(beta.^Tlag).*(1-pi).*(1-beta))./(pi.*(beta.^Tlag).*(1-beta));

deforest = zeros(1,Tmax);
emissions = zeros(1,Tmax);
DiffProfit = zeros(N,N);  % This will store the differences in Agribusinessess
                          % function psi in the model
figure;                % Starting animation
h = imagesc(D(:,:,1)); % Create heatmap
colorbar;
for t = 1:Tmax-1
    % ===== Step 1: calculating agribusiness diff profits =====
    B = (Agribusiness(:,:,t) == 1);                % current agribusiness map (logical)
    
    % -- find 8-connected clusters (requires Image Processing Toolbox) --
    CC = bwconncomp(B, 8);
    clusters = CC.PixelIdxList;                    % cell array of linear indices
    
    % -- border neighbors of each cluster via convolution (4-neighborhood) --
    K8 = ones(3); K8(2,2)=0;                    % Von Neumann neighborhood
    neighbors = cell(size(clusters));
    for i = 1:numel(clusters)
        M = false(N,N); M(clusters{i}) = true;
        borderMask = conv2(double(M), K8, 'same') > 0 & ~M;   % adjacent but not in cluster
        neighbors{i} = find(borderMask);                      % linear indices
    end
    
    % -- compute DiffProfit only on border cells not already agribusiness --
    DiffProfit = zeros(N,N);
    for i = 1:numel(neighbors)
        clusterIdx = clusters{i};
        ContAgri   = numel(clusterIdx);                       % size of cluster
    
        % profit() should accept linear indices; term1/term2 are scalars
        term1 = sum(profit(A, ContAgri+1, clusterIdx));
        term2 = sum(profit(A, ContAgri,   clusterIdx));
    
        nb = neighbors{i};
        nb = nb(~B(nb));                                      % exclude cells already in agribusiness
    
        % fill marginal diff profit on border candidates
        DiffProfit(nb) = term1 - term2 + profit(A, ContAgri+1, nb);
    end
    
    % ===== Step 2 (NxN): farmers' profit and price R =====
    q1
    FarmersProfit = max(q1.*F + q2.*wage - q3.*A - wage/(1-beta), (A - wage)/(1-beta));
    R(:,:,t) = w1.*DiffProfit + (1-w1).*FarmersProfit;              % adjust weights as desired
    
    % ===== Farmers' decisions (NxN) =====
    Vu = value_u(A, R(:,:,t), F);
    Vp = value_p(A, F);
    
    Du = (L == 0) & (D(:,:,t) ~= 1) & (Agribusiness(:,:,t) ~= 1) & (Vu > 0);
    Dp = (L == 1) & (D(:,:,t) ~= 1) & (Agribusiness(:,:,t) ~= 1) & (Vp > 0);
    
    % ===== State updates (NxN) =====
    D(:,:,t+1) = D(:,:,t) + Du + Dp;
    O(:,:,t+1) = O(:,:,t) + D(:,:,t) .* O_aux;               % O_aux must be N x N
    
    Candidates = D(:,:,t) - Agribusiness(:,:,t);              % plots that can gain rights
    PP = (Candidates == 1) & (rand(N,N) < pi) & (L == 0);      % earn property rights
    NewAgri = (PP == 1) & (DiffProfit >= R(:,:,t));
    Agribusiness(:,:,t+1) = Agribusiness(:,:,t) + NewAgri;
    
    % (optional) visualization + aggregates
    set(h, 'CData', D(:,:,t));                                % heatmap update
    pause(0.1); % Pause for smooth animation
    deforest(t)  = sum(D(:,:,t), 'all');
    emissions(t) = sum(O(:,:,t), 'all');
end

% === Final Heatmaps ===
figure;
subplot(1,2,1)

% --- Base map (D) ---
imagesc(D(:,:,1));
colormap([0 1 0; 1 0 0]);   % green for 0, red for 1
caxis([0 1]);
hold on;

% --- Overlay (L) only where D==0 ---
overlayRGB = zeros([size(L), 3]);  % initialize RGB image
overlayRGB(:,:,3) = 1;             % blue channel = 1

% Alpha mask: 1 only where L==1 and D==0
alphaMask = double(L == 1 & squeeze(D(:,:,1)) == 0);

% Expand alpha mask to match image size if needed
hOverlay = imagesc(overlayRGB);
set(hOverlay, 'AlphaData', alphaMask);

hold off;
colorbar;
title('Initial Deforestation');


subplot(1,2,2)
imagesc(D(:,:,end)); colormap([0 1 0; 1 0 0]);   % green for 0, red for 1
clim([0 1]);
hold on;

% --- Overlay (L) only where D==0 ---
overlayRGB = zeros([size(L), 3]);  % initialize RGB image
overlayRGB(:,:,3) = 1;             % blue channel = 1

% Alpha mask: 1 only where L==1 and D==0
alphaMask = double(L == 1 & squeeze(D(:,:,end)) == 0);

% Expand alpha mask to match image size if needed
hOverlay = imagesc(overlayRGB);
set(hOverlay, 'AlphaData', alphaMask);

hold off;
colorbar;
title('Initial Deforestation with L overlay');

figure;
subplot(1,2,1)
plot(1:Tmax-1, deforest(1:Tmax-1), '-o', 'LineWidth', 1.5);
xlabel('Time');
ylabel('Accumulated Deforestation');
title('Deforestation over Time');
grid on;

subplot(1,2,2)
plot(1:Tmax-1, emissions(1:Tmax-1), '-o', 'LineWidth', 1.5);
xlabel('Time');
ylabel('Emissions');
title('Emissions over Time');
grid on;