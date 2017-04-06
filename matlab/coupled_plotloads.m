function coupled_plotloads(mesh, loads)

n = size(mesh,1)/2;  % number of spanwise points

% Rearrange mesh array to order the wing points for drawing
tmpmesh = zeros(size(mesh,1)+1,size(mesh,2));
tmpmesh(1:n,:) = mesh(1:n,:);
tmpmesh(n+1:end-1,:) = flipud(mesh(n+1:end,:));
tmpmesh(end,:) = mesh(1,:);
% mesh = tmp;
M1 = tmpmesh(:,1);  M2 = tmpmesh(:,2);  M3 = tmpmesh(:,3);

fem_origin = 0.35;   % fem origin at 0.35*chord

% mesh(end,:) = [];  % remove first point at end
Lx = zeros(n,1); Ly = zeros(n,1); Lz = zeros(n,1);
for i = 1:n  % x cooridates for force,moments
    Lx(i) = fem_origin*(mesh(i+n,1) - mesh(i,1)) + mesh(i,1); 
end

% interpolate to get y and z coordinates
for i = 1:n
    x = [mesh(i,1),mesh(i+n,1)];
    y = [mesh(i,2),mesh(i+n,2)];
    z = [mesh(i,3),mesh(i+n,3)];
    Ly(i) = interp1(x,y,Lx(i));
    Lz(i) = interp1(x,z,Lx(i));
end

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create mesh plot
mp = plot3(M1,M2,M3,'-bo','LineWidth',1.5,'MarkerFaceColor','b');

% Create labels
xlabel(axes1,'X-axis (chord) [m]');
ylabel(axes1,'Y-axis (span) [m]');
zlabel(axes1,'Z-axis (elevation) [m]');
ylim(axes1,[-40,40]);
fill3(M1,M2,M3,[204 229 255]./255);  % fill in top of wing
fill3(M1,M2,M3-.001,[0 76 153]./255);  % fill in underside of wing
plot3(Lx,Ly,Lz,'-ro','LineWidth',1.5,'MarkerFaceColor','r','MarkerSize',5); % vector origins
% scF=5e-6; scM=1e-5;
scF = 1; scM = 1;
Fu = loads(:,1)*scF; Fv = loads(:,2)*scF; Fw = loads(:,3)*scF;
Mu = loads(:,4)*scM; Mv = loads(:,5)*scM; Mw = loads(:,6)*scM;

qf = quiver3(Lx,Ly,Lz,Fu,Fv,Fw,0.5,'Linewidth',2,'Color','r'); % Forces
qm = quiver3(Lx,Ly,Lz,Mu,Mv,Mw,2,'Linewidth',2,'Color','m');   % Moments
ql = legend([qf qm mp],'Load Forces','Load Moments','Mesh Points');
ql.FontSize = 12;

view(axes1,[78.1 10.8]);
box(axes1,'on');
grid(axes1,'on');
