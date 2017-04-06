function coupled_plotmesh(mesh,origmesh)

n = size(mesh,1)/2;  % number of spanwise points

% Rearrange mesh array to order the wing points for drawing
tmpmesh = zeros(size(mesh,1)+1,size(mesh,2));
tmpmesh(1:n,:) = mesh(1:n,:);
tmpmesh(n+1:end-1,:) = flipud(mesh(n+1:end,:));
tmpmesh(end,:) = mesh(1,:);
M1 = tmpmesh(:,1);  M2 = tmpmesh(:,2);  M3 = tmpmesh(:,3);

% Rearrange original mesh array to order the wing points for drawing
tmpmesh = zeros(size(origmesh,1)+1,size(origmesh,2));
tmpmesh(1:n,:) = origmesh(1:n,:);
tmpmesh(n+1:end-1,:) = flipud(origmesh(n+1:end,:));
tmpmesh(end,:) = origmesh(1,:);
OM1 = tmpmesh(:,1);  OM2 = tmpmesh(:,2);  OM3 = tmpmesh(:,3);

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create mesh plot
mp = plot3(M1,M2,M3,'-bo','LineWidth',1.5,'MarkerFaceColor','b');
op = plot3(OM1,OM2,OM3,'-ro','LineWidth',1.5,'MarkerFaceColor','r');

% Create labels
xlabel(axes1,'X-axis (chord) [m]');
ylabel(axes1,'Y-axis (span) [m]');
zlabel(axes1,'Z-axis (elevation) [m]');
zlimits = .075;
zlim(axes1,[-zlimits*.5,zlimits]);
ylim(axes1,[-40,40]);
dw = fill3(M1,M2,M3,[204 229 255]./255);  % fill in wing light blue
% fill3(M1,M2,M3-.0005,[0 76 153]./255);  % fill in underside of wing
ow = fill3(OM1,OM2,OM3,[255 102 102]./255);  % fill in orig wing light red
ml = legend([mp op],'Displaced Mesh','Initial Mesh');
ml.FontSize = 12;

view(axes1,[78.1 10.8]);
box(axes1,'on');
grid(axes1,'on');
