function h = plot3dBalloon(R, AZ, EL, center, varargin)
%% 1. Parsing e validazione degli argomenti in ingresso
p = inputParser;
addRequired(p, 'R', @isnumeric);
addRequired(p, 'AZ', @isnumeric);
addRequired(p, 'EL', @isnumeric);
addRequired(p, 'center', @(x) isnumeric(x) && numel(x) == 3);
addParameter(p, 'Units', 'deg', @(x) ischar(x) || isstring(x));
addParameter(p, 'FaceAlpha', 0.85, @(x) isnumeric(x) && x>=0 && x<=1);
addParameter(p, 'EdgeColor', 'none');
addParameter(p, 'Parent', [], @(x) isempty(x) || isgraphics(x, 'axes'));
% nuovi parametri per i limiti degli assi (1x2 numeric arrays)
addParameter(p, 'XLim', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
addParameter(p, 'YLim', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
addParameter(p, 'ZLim', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
parse(p, R, AZ, EL, center, varargin{:});
units = lower(p.Results.Units);
alphaVal = p.Results.FaceAlpha;
edgeCol = p.Results.EdgeColor;
ax = p.Results.Parent;
xlimUser = p.Results.XLim;
ylimUser = p.Results.YLim;
zlimUser = p.Results.ZLim;
if isempty(ax)
    ax = gca;
end

%% 2. Gestione delle griglie angolari (1D -> 2D meshgrid)
% Se AZ ed EL sono vettori 1D, crea le matrici 2D corrispondenti
if isvector(AZ) && isvector(EL) && ~isvector(R)
    [AZ_grid, EL_grid] = meshgrid(AZ, EL);
else
    AZ_grid = AZ;
    EL_grid = EL;
end

% Conversione in radianti se gli angoli sono espressi in gradi
if strcmp(units, 'deg')
    AZ_rad = deg2rad(AZ_grid);
    EL_rad = deg2rad(EL_grid);
else
    AZ_rad = AZ_grid;
    EL_rad = EL_grid;
end

%% 3. Trasformazione da coordinate Sferiche a Cartesiane
% MATLAB sph2cart usa: Azimuth (rad), Elevation (rad da -pi/2 a pi/2), Raggio
[X_local, Y_local, Z_local] = sph2cart(AZ_rad, EL_rad, R);

% Traslazione del balloon nel centro di espansione desiderato [x0, y0, z0]
X = X_local + center(1);
Y = Y_local + center(2);
Z = Z_local + center(3);

%% 4. Rendering grafico del Balloon
% Mantiene lo stato di hold precedente per consentire sovrapposizioni
wasHolding = ishold(ax);
hold(ax, 'on');
% Tracciamento della superficie 3D (la mappa colori rappresenta R)
h = surf(ax, X, Y, Z, R, ...
    'FaceAlpha', alphaVal, ...
    'EdgeColor', edgeCol, ...
    'FaceColor', 'interp');

% Configurazione dello scenario 3D se è il primo elemento tracciato
if ~wasHolding
    axis(ax, 'equal');
    grid(ax, 'on');
    xlabel(ax, 'X [m]');
    ylabel(ax, 'Y [m]');
    zlabel(ax, 'Z [m]');
    view(ax, 3);
    colormap(ax, 'jet');
    colorbar(ax);
    lighting(ax, 'gouraud');
    camlight(ax, 'headlight');
end

% Applica i limiti degli assi se forniti dall'utente
% Nota: impostiamo i limiti DOPO axis equal per forzarli
if ~isempty(xlimUser)
    xlim(ax, xlimUser);
end
if ~isempty(ylimUser)
    ylim(ax, ylimUser);
end
if ~isempty(zlimUser)
    zlim(ax, zlimUser);
end

% ripristina lo stato di hold precedente
if ~wasHolding
    hold(ax, 'off');
end
end