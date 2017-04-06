% Draw the CRM wing used in OpenAeroStruct

% points
p = zeros(6,3);

p(1,:) = [904.294, 0, 174.126]; % x,y,z of LE
p(6,:) = p(1,:) + [536.181, 0, 0]; % choard length

p(2,:) = [1208.203, 404.864, 180.359];
p(5,:) = p(2,:) + [299.317, 0, 0];

p(3,:) = [1780.737, 1156.753, 263.827];
p(4,:) = p(3,:) + [107.4, 0, 0];


% Zero out the z-components
p(:,3) = 0;

p

% Draw diagram
wing = fill(p(:,1),p(:,2),'white');
% rotate(wing, [1, 0, 0], 90)
