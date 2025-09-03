function hEst = plot_estados_top(ax, estados, lnclr)
    hEst = gobjects(0);
    hold(ax,'on');
    % Prefer mapshow when available — it handles many shapes automatically
    if exist('mapshow','file') == 2
        try
            hEst = mapshow(ax, estados, 'DisplayType','line', 'Color', lnclr, 'LineWidth',1.0);
            uistack(hEst,'top'); return;
        catch
            % fall through to manual cases
        end
    end
    % Manual cases
    if isstruct(estados)
        if isfield(estados,'X') && isfield(estados,'Y')
            for k = 1:numel(estados)
                hEst(end+1) = plot(ax, estados(k).X, estados(k).Y, '-', 'Color', lnclr, 'LineWidth',1.0); %#ok<AGROW>
            end
        elseif isfield(estados,'Lon') && isfield(estados,'Lat')
            for k = 1:numel(estados)
                hEst(end+1) = plot(ax, estados(k).Lon, estados(k).Lat, '-', 'Color', lnclr, 'LineWidth',1.0); %#ok<AGROW>
            end
        end
    elseif istable(estados) && any(strcmpi('Shape', estados.Properties.VariableNames))
        S = estados.Shape;
        % table of shapes -> plot each
        for k = 1:height(estados)
            try
                hEst(end+1) = plot(ax, S(k), 'EdgeColor', lnclr, 'FaceColor', 'none', 'LineWidth',1.0); %#ok<AGROW>
            catch
                % If Shape stores polyshape in cells:
                try
                    shp = S{k};
                    hEst(end+1) = plot(ax, shp, 'EdgeColor', lnclr, 'FaceColor','none','LineWidth',1.0); %#ok<AGROW>
                catch
                end
            end
        end
    elseif isa(estados,'polyshape')
        hEst = plot(ax, estados, 'FaceColor','none', 'EdgeColor', lnclr, 'LineWidth',1.0);
    else
        % Last-ditch attempt (numeric Nx2, or graphics-ready)
        try
            hEst = plot(ax, estados, 'Color', lnclr, 'LineWidth',1.0);
        catch
            % If you still hit this, uncomment next line to inspect the type:
            % disp(class(estados));
        end
    end
    if ~isempty(hEst), uistack(hEst,'top'); end
end