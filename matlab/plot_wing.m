function [] = plot_wing(wing_twist_cp, wing_thickness_cp, wing_taper, wing_chord_cp)

% import OpenAeroStruct python module
OAS_PATH = py.os.path.abspath('../..');
P = py.sys.path;
if count(P,OAS_PATH) == 0
    insert(P,int32(0),OAS_PATH);
end
            
% get mesh, twist, thickness values for design variables
prob_dict = struct;
prob_dict.type = 'aerostruct';
prob_dict.with_viscous = true;
prob_dict.optimize = false;
prob_dict.record_db = false;  % using sqlitedict locks a process
prob_dict.print_level = 0;
prob_dict.alpha = 0.;
OAS_prob = py.OpenAeroStruct.run_classes.OASProblem(prob_dict);
% Create a dictionary to store options about the surface
surf_dict = struct;
surf_dict.name = 'wing';
surf_dict.num_y = 7;
surf_dict.num_x = 2;
surf_dict.wing_type = 'CRM';
surf_dict.CD0 = 0.015;
surf_dict.symmetry = true;
surf_dict.num_twist_cp = 2;
surf_dict.num_thickness_cp = 2;
surf_dict.num_chord_cp = 1;
surf_dict.exact_failure_constraint = true;
surf_dict.span_cos_spacing = 0.5;
OAS_prob.add_surface(surf_dict);
% Multiple lifting surfaces
surf_dict = struct;
surf_dict.name = 'tail';
surf_dict.num_y = 7;
surf_dict.num_x = 2;
surf_dict.span = 20.;
surf_dict.root_chord = 5.;
surf_dict.wing_type = 'rect';
surf_dict.offset = [50., 0., 5.];
surf_dict.twist_cp = -9.5;
surf_dict.exact_failure_constraint = true;
OAS_prob.add_surface(surf_dict)
OAS_prob.add_desvar('wing.twist_cp');
OAS_prob.add_desvar('wing.thickness_cp');
OAS_prob.add_desvar('wing.taper');
OAS_prob.add_desvar('wing.chord_cp');
OAS_prob.setup()
% Actually run the problem
input = {'wing.twist_cp',wing_twist_cp,'wing.thickness_cp',wing_thickness_cp,...
    'wing.taper',wing_taper,'wing.chord_cp',wing_chord_cp,'matlab',true};
data = struct(OAS_prob.run(pyargs(input{:})));
wingmesh = np2mat(data.wing_mesh);
thickness = np2mat(data.wing_thickness);
twist = np2mat(data.wing_twist);
vonmises = np2mat(data.wing_vonmises);
vonmises = max(vonmises,[],2);
yield = struct(OAS_prob.surfaces{1}).yield;

% initial CRM wing
wingmesh_init_x = ...
    [4.5230719800000003e+01, 4.0003449900510176e+01, 3.2179945800000006e+01, 2.2969067599999999e+01;
    4.7958679800000006e+01, 4.4411817603600596e+01, 3.9102089266666667e+01, 3.6588064999999993e+01];
    
wingmesh_init_y = ...
    [-2.9381526199999996e+01, -2.2516417064784868e+01, -1.2242300466666668e+01, -8.9954746742019987e-16;
    -2.9381526199999996e+01, -2.2516417064784868e+01, -1.2242300466666668e+01, -8.9954746742019987e-16];
    
wingmesh_init_z = ...
    [6.7012057999999994e+00, 5.5812274808896820e+00, 4.6522809333333335e+00, 4.4228003999999999e+00;
    6.7012057999999994e+00, 5.5812274808896820e+00, 4.6522809333333335e+00, 4.4228003999999999e+00];

wingmesh_init = cat(3,wingmesh_init_x, wingmesh_init_y, wingmesh_init_z);

% [[[ 4.5230719800000003e+01 -2.9381526199999996e+01
%     6.7012057999999994e+00]
%   [ 4.0003449900510176e+01 -2.2516417064784868e+01
%     5.5812274808896820e+00]
%   [ 3.2179945800000006e+01 -1.2242300466666668e+01
%     4.6522809333333335e+00]
%   [ 2.2969067599999999e+01 -8.9954746742019987e-16
%     4.4228003999999999e+00]]
% 
%  [[ 4.7958679800000006e+01 -2.9381526199999996e+01
%     6.7012057999999994e+00]
%   [ 4.4411817603600596e+01 -2.2516417064784868e+01
%     5.5812274808896820e+00]
%   [ 3.9102089266666667e+01 -1.2242300466666668e+01
%     4.6522809333333335e+00]
%   [ 3.6588064999999993e+01 -8.9954746742019987e-16
%     4.4228003999999999e+00]]]

% mesh matrix:
% mesh(nx, ny, [x,y,z])
%    mesh=(LE -> TE, Mid -> Tip, [x,y,z])

x = wingmesh(:,:,1);
y = -wingmesh(:,:,2);
z = wingmesh(:,:,3);

nx = size(wingmesh,1);  % number of chordwise points
ny = size(wingmesh,2);  % number of spanwise points
lw = 2;

% span = m_vals[0, -1, 1] - m_vals[0, 0, 1]
% rel_span = (m_vals[0, :, 1] - m_vals[0, 0, 1]) * 2 / span - 1
% span_diff = ((m_vals[0, :-1, 1] + m_vals[0, 1:, 1]) / 2 - m_vals[0, 0, 1]) * 2 / span - 1
span = wingmesh(1,end,2) - wingmesh(1,1,2);
rel_span = (wingmesh(1,:,2)-wingmesh(1,1,2))/span;
span_diff = -(wingmesh(1,1:end-1,2) + wingmesh(1,2:end,2))/(2*span);


fig = figure(1);
ax = gca;
hold(ax,'on');
% plot optimized wing
for i = [1,nx]
    plot(y(i,:),x(i,:),'-k','LineWidth',lw)
end
p1 = plot(y(:,[1,end]),x(:,[1,end]),'-k','LineWidth',lw);
for i = [2,ny-1]
    plot(y(:,i),x(:,i),'-k','LineWidth',lw)
end
% plot initial wing
yinit = -wingmesh_init_y;
xinit = wingmesh_init_x;
for i = [1,nx]
    plot(yinit(i,:),xinit(i,:),'--k','LineWidth',lw/2)
end
p2 = plot(yinit(:,[1,end]),xinit(:,[1,end]),'--k','LineWidth',lw/2);
for i = [2,ny-1]
    plot(yinit(:,i),xinit(:,i),'--k','LineWidth',lw/2)
end
legend([p1(1),p2(1)],{'Optimized','Initial'});
ax.YAxis.Direction = 'reverse';
ax.YAxis.Label.String = 'Chord [m]';
ax.XAxis.Label.String = 'Span [m]';
ax.XAxis.Limits = [-2.5,max(y(:,1))+2.5];

% xlabel('Y');
% ylabel('X');

% Plot the spar thickness overlay
% lw_thick = 5;
% n = 20; % number of interpolated colors to draw
% % get individual colors from parula colormap
% map = colormap('parula');
% lb_color = 0.0001;
% ub_color = 0.06;
% caxis([lb_color, ub_color]);
% y_thick = (y(1,2:end)+y(1,1:end-1))/2;  % y-coord for thickness values
% y_norm = y_thick/max(y_thick);
% y_q = linspace(min(y_thick),max(y_thick),n);
% thick_q = interp1(y_thick,thickness,y_q);
% % add boundary points to thickness values with extrapolation
% thick_q_norm = (thick_q-min(thick_q))/(max(thick_q)-min(thick_q)); % normalized thickness
% % thick_q_norm = (thick_q-lb_color)./(ub_color-lb_color); % normalized thickness
% thick_q2 = thick_q_norm*(size(map,1)-1) + 1;
% % y_norm = (y(1,1:end-1)-y(1,1))/(y(1,end)-y(1,1));
% % xq = linspace(thickness(1),thickness(end),n);
% idx_norm = (1:size(map,1))/(size(map,1));
% c_q = interp1(1:size(map,1),map,thick_q2);
% 
% x_loc = x(1, :)*.5 + x(end,:)*.5;
% x_q = interp1(y(1,:),x_loc,y_q);
% % y_spar = 
% for i = 1:length(thick_q)-1
%     plot(y_q(i:i+1),x_q(i:i+1),'-','Color',c_q(i,:),'LineWidth',lw_thick);
% end
% hold off
% for i = 1:length(x_loc)-1
%     spar = plot3(x_loc(i:i+2), y(0,i:i+2),...
%         c=parula.parula_map(thickness[i]), lw=5, solid_capstyle='butt')
% end

% Plot thickness
fig2 = figure(2);
ax = gca;
plot(span_diff, thickness,'LineWidth',lw)
ax.YAxis.Label.String = 'Spar Thickness [m]';
ax.XAxis.Label.String = 'Relative Span';
ax.XAxis.Limits = [0,1];

% Plot twist
fig2 = figure(3);
ax = gca;
plot(rel_span, twist,'LineWidth',lw)
ax.YAxis.Label.String = 'Twist [deg]';
ax.XAxis.Label.String = 'Relative Span';
ax.XAxis.Limits = [0,1];

% Plot vonmises
fig2 = figure(4);
ax = gca;
hold(ax,'on');
plot(span_diff, vonmises,'-b','LineWidth',lw)
plot([0,1],[yield, yield],'--r','LineWidth',lw/1.5);
ax.YAxis.Label.String = 'Von Mises';
ax.XAxis.Label.String = 'Relative Span';
ax.XAxis.Limits = [0,1];
ax.YAxis.Limits = [min(vonmises)*0.95, max(vonmises)*1.05];
legend({'Von Mises Stress','Yield Stress'});


% % Rearrange mesh array to order the wing points for drawing
% tmpmesh = zeros(size(mesh,1)+1,size(mesh,2));
% tmpmesh(1:n,:) = mesh(1:n,:);
% tmpmesh(n+1:end-1,:) = flipud(mesh(n+1:end,:));
% tmpmesh(end,:) = mesh(1,:);
% M1 = tmpmesh(:,1);  M2 = tmpmesh(:,2);  M3 = tmpmesh(:,3);
% 
% % Rearrange original mesh array to order the wing points for drawing
% tmpmesh = zeros(size(origmesh,1)+1,size(origmesh,2));
% tmpmesh(1:n,:) = origmesh(1:n,:);
% tmpmesh(n+1:end-1,:) = flipud(origmesh(n+1:end,:));
% tmpmesh(end,:) = origmesh(1,:);
% OM1 = tmpmesh(:,1);  OM2 = tmpmesh(:,2);  OM3 = tmpmesh(:,3);
% 
% % Create figure
% figure1 = figure;
% 
% % Create axes
% axes1 = axes('Parent',figure1);
% hold(axes1,'on');
% 
% % Create mesh plot
% mp = plot3(M1,M2,M3,'-bo','LineWidth',1.5,'MarkerFaceColor','b');
% op = plot3(OM1,OM2,OM3,'-ro','LineWidth',1.5,'MarkerFaceColor','r');
% 
% % Create labels
% xlabel(axes1,'X-axis (chord) [m]');
% ylabel(axes1,'Y-axis (span) [m]');
% zlabel(axes1,'Z-axis (elevation) [m]');
% zlimits = .075;
% zlim(axes1,[-zlimits*.5,zlimits]);
% ylim(axes1,[-40,40]);
% dw = fill3(M1,M2,M3,[204 229 255]./255);  % fill in wing light blue
% % fill3(M1,M2,M3-.0005,[0 76 153]./255);  % fill in underside of wing
% ow = fill3(OM1,OM2,OM3,[255 102 102]./255);  % fill in orig wing light red
% ml = legend([mp op],'Displaced Mesh','Initial Mesh');
% ml.FontSize = 12;
% 
% view(axes1,[78.1 10.8]);
% box(axes1,'on');
% grid(axes1,'on');
end
