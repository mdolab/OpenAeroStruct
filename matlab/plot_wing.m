data = output;
if isstruct(data) 
    if isa(data.wing_mesh,'py.numpy.ndarray')
        mesh = np2mat(data.wing_mesh);
    end
    if isa(data.wing_thickness, 'py.numpy.ndarray')
        thickness = np2mat(data.wing_thickness);
    end
else
end

% mesh matrix:
% mesh(nx, ny, [x,y,z])
%    mesh=(LE -> TE, Mid -> Tip, [x,y,z])

x = mesh(:,:,1);
y = -mesh(:,:,2);
% z = mesh(:,:,3);

nx = size(mesh,1);  % number of chordwise points
ny = size(mesh,2);  % number of spanwise points
lw = 2;
hold on
for i = [1,nx]
    plot(y(i,:),x(i,:),'-k','LineWidth',lw)
end
plot(y(:,[1,end]),x(:,[1,end]),'-k','LineWidth',lw)
ax = gca;
ax.YAxis.Direction = 'reverse';
ax.YAxis.Label.String = 'X';
ax.XAxis.Label.String = 'Y';
ax.XAxis.Limits = [-2.5,max(y(:,1))+2.5];
% xlabel('Y');
% ylabel('X');

% Plot the spar thickness overlay
n = 100; % number of interpolated colors to draw
% get individual colors from parula colormap
map = colormap('parula');
lb_color = 0.0001;
ub_color = 0.06;
caxis([lb_color, ub_color]);
y_thick = (y(1,2:end)+y(1,1:end-1))/2;  % y-coord for thickness values
y_norm = y_thick/max(y_thick);
y_q = linspace(min(y_thick),max(y_thick),n);
thick_q = interp1(y_thick,thickness,y_q);
% add boundary points to thickness values with extrapolation
thick_q_norm = (thick_q-min(thick_q))/(max(thick_q)-min(thick_q)); % normalized thickness
% thick_q_norm = (thick_q-lb_color)./(ub_color-lb_color); % normalized thickness
thick_q2 = thick_q_norm*(size(map,1)-1) + 1;
% y_norm = (y(1,1:end-1)-y(1,1))/(y(1,end)-y(1,1));
% xq = linspace(thickness(1),thickness(end),n);

idx_norm = (1:size(map,1))/(size(map,1));
c_q = interp1(1:size(map,1),map,thick_q2);

x_loc = x(1, :)*.5 + x(end,:)*.5;
x_q = interp1(y(1,:),x_loc,y_q);
% y_spar = 
for i = 1:length(thick_q)-1
    plot(y_q(i:i+1),x_q(i:i+1),'-','Color',c_q(i,:),'LineWidth',10);
end
hold off
% for i = 1:length(x_loc)-1
%     spar = plot3(x_loc(i:i+2), y(0,i:i+2),...
%         c=parula.parula_map(thickness[i]), lw=5, solid_capstyle='butt')
% end
return
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
