function SolveDeforestation1D(params,A,F,O,D,L,Agribusiness,O_aux,w1,animation)

q1 = (1-params.beta.*(1-params.pi))./(params.pi.*(params.beta.^params.Tlag));
q2 = (1-params.beta.*(1-params.pi).*(1+params.pi.*(params.beta.^params.Tlag)))./(params.pi.*(params.beta.^params.Tlag).*(1-params.beta));
q3 = (params.alpha.*(1-params.beta.^params.Tlag)+(params.beta.^params.Tlag).*(1-params.pi).*(1-params.beta))./(params.pi.*(params.beta.^params.Tlag).*(1-params.beta));

deforest = zeros(1,params.Tmax);
emissions = zeros(1,params.Tmax);
DiffProfit = zeros(1,params.N);  % This will store the differences in Agribusinessess

if animation == true% function psi in the model
    figure;                % Starting animation
    h = imagesc(D(1,:,1)); % Create heatmap
    colorbar;
    set(h, 'CData', D(1,:,1)); % Update heatmap data
    pause(1.5); % Pause for smooth animation
end
for t = 1:params.Tmax-1
    % Calculating land prices R
    % Step 1: calcuting agribusiness diff profits
    clusters = {};
    isOne = (Agribusiness(1,:,t) == 1);
    d = diff([0 isOne 0]);
    startIdx = find(d == 1);
    endIdx   = find(d == -1) - 1;
    
    for k = 1:length(startIdx)
        clusters{k} = startIdx(k):endIdx(k);
    end
    
    left  = cellfun(@(c) min(c(:)), clusters);
    right = cellfun(@(c) max(c(:)), clusters);
    
    cand  = [left-1; right+1];              
    valid = cand >= 1 & cand <= params.N;          % in-bounds mask
    
    neighbors = arrayfun(@(i) cand(valid(:,i), i).', ...
                         1:numel(left), 'UniformOutput', false);
    
     for i =1:numel(neighbors)
         cluster = clusters{i};
         ContAgri = numel(cluster); % Counting the S
         term1 = sum(profit(params,A,ContAgri+1,cluster));
         term2 = sum(profit(params,A,ContAgri,cluster));
         for k = 1:numel(neighbors{i})
                if Agribusiness(1,neighbors{i}(k),t) == 1   % we just need to compute this for plots    
                                                            % which have not been sold to a Agri
                           continue;
                end
        DiffProfit(1,neighbors{i}(k)) = term1 - term2 + profit(params,A,ContAgri+1,neighbors{i}(k));   
         end
     end
        % Step 2: calculate the farmers profit when selling land
     FarmersProfit = max(q1*F+q2*params.wage-q3*A-params.wage/(1-params.beta), (A-params.wage)/(1-params.beta));

    R(1,:,t) = w1.*DiffProfit + (1-w1).*FarmersProfit;
     
    % Farmers
    Vu = value_u(params,A,R(1,:,t),F);
    Vp = value_p(params,A,F);

    Du = L == 0 & D(1,:,t) ~= 1 & Agribusiness(1,:,t) ~= 1 & Vu > 0; 
    Dp = L == 1 & D(1,:,t) ~= 1 & Agribusiness(1,:,t) ~= 1 & Vp > 0;
    
    % Updating Deforestation
    D(1,:,t+1) = D(1,:,t) + Du+Dp;
    % Updating Emissions
    O(1,:,t+1) = O(1,:,t) +  D(1,:,t) .* O_aux;
    % Updating the Agribusiness
    Candidates = D(1,:,t) - Agribusiness(1,:,t); %plots that can gain property rights
    PP = Candidates ==1 & rand(1,params.N) < params.pi & L==0; %plots that have earned property rights
    NewAgri = PP==1 & DiffProfit >= R(1,:,t);
    Agribusiness(1,:,t+1) = Agribusiness(1,:,t) + NewAgri;
    if animation == true
        set(h, 'CData', D(1,:,t)); % Update heatmap data
        pause(0.1); % Pause for smooth animation
    end
    deforest(t) = sum(D(1,:,t));
    emissions(t) = sum(O(1,:,t));
end

% === Final Heatmaps ===
% Figure 1: Initial x Final Feforestiation
figure;
subplot(1,2,1)

% --- Base map (D) ---
imagesc(D(1,:,1));
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
alphaMask = double(L == 1 & squeeze(D(1,:,1)) == 0);

% Expand alpha mask to match image size if needed
hOverlay = imagesc(overlayRGB);
set(hOverlay, 'AlphaData', alphaMask);

hold off;
title('Initial Deforestation');


subplot(1,2,2)
imagesc(D(1,:,end)); colormap([0 1 0; 1 0 0]);   % green for 0, red for 1
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
alphaMask = double(L == 1 & squeeze(D(1,:,end)) == 0);

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