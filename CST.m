clc,
clear,

format long g

%% INPUTS %%

nodes = readmatrix("nodes.txt");
nodes = nodes(:,2:end);
nnd = length(nodes); % Number of nodes %

elements = readmatrix("elements.txt");
elements = elements(:,2:end);
nel = length(elements); % Number of elements %

bc = readmatrix("boundaries.txt");

nodal_loads = readmatrix("nloads.txt");
tractions = readmatrix("tractions.txt");
vol_loads = readmatrix("vloads");

%% MATERIAL PROPERTIES %%

E = 70e6;
nu = 0.25;
thickness = 0.02;
plane = 0; % 0 if plane stress , 1 if plane strain

if plane==0
    D = (E/(1-nu^2)) * [1 nu 0;
                        nu 1 0;
                         0 0 (1-nu)/2];
elseif plane==1
    D = (E/(1+nu)/(1-2*nu)) * [1-nu -nu 0;
                               -nu 1-nu 0;
                               0 0 (1-2*nu)/2];
end

%% ASSEMBLY OF THE FORCES VECTOR %%

F = zeros(2*nnd,1); % Global force vector

for i=1:size(nodal_loads,1)

    Fnode=nodal_loads(i,1);
    Fmag=nodal_loads(i,2);
    Fangle=nodal_loads(i,3)/180*pi;

    fx=Fmag*cos(Fangle);
    fy=Fmag*sin(Fangle);

    F(2*Fnode-1,1) = F(2*Fnode-1) + fx;
    F(2*Fnode,1) = F(2*Fnode) + fy;
end

for i=1:size(tractions,1)

    first_node=tractions(i,1);
    end_node=tractions(i,2);
    traction_mag=tractions(i,3);
    traction_angle=tractions(i,4)/180*pi;

    x1=nodes(first_node,1);
    x2=nodes(end_node,1);
    Dx=x2-x1;

    y1=nodes(first_node,2);
    y2=nodes(end_node,2);
    Dy=y2-y1;

    L=sqrt(Dx^2+Dy^2);

    Tx=0.5 * thickness * traction_mag * L * cos(traction_angle);
    Ty=0.5 * thickness * traction_mag * L * sin(traction_angle);

    F(2*first_node-1,1) = F(2*first_node-1) + Tx;
    F(2*first_node,1) = F(2*first_node) + Ty;
    F(2*end_node-1,1) = F(2*end_node-1) + Tx;
    F(2*end_node,1) = F(2*end_node) + Ty; 
end

for i=1:size(vol_loads,1)

    vol_element=vol_loads(i,1);
    vol_mag=vol_loads(i,2);
    vol_angle=vol_loads(i,3)/180*pi;

    n1=elements(vol_element,1);
    n2=elements(vol_element,2);
    n3=elements(vol_element,3);

    x1=nodes(n1,1);
    x2=nodes(n2,1);
    x3=nodes(n3,1);

    y1=nodes(n1,2);
    y2=nodes(n2,2);
    y3=nodes(n3,2);

    A=(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2))/2;

    bx=1/3 * vol_mag * thickness * A * cos(vol_angle);
    by=1/3 * vol_mag * thickness * A * sin(vol_angle);

    F(2*n1-1,1) = F(2*n1-1) + bx;
    F(2*n1,1) = F(2*n1) + by;
    F(2*n2-1,1) = F(2*n2-1) + bx;
    F(2*n2,1) = F(2*n2) + by;
    F(2*n3-1,1) = F(2*n3-1) + bx;
    F(2*n3,1) = F(2*n3) + by;
end

%% ASSEMBLY OF THE GLOBAL STIFFNESS MATRIX %%

K = zeros(2*nnd,2*nnd); % Global stiffness matrix

for i=1:size(elements,1)

    n1=elements(i,1);
    n2=elements(i,2);
    n3=elements(i,3);

    x1=nodes(n1,1);
    x2=nodes(n2,1);
    x3=nodes(n3,1);

    y1=nodes(n1,2);
    y2=nodes(n2,2);
    y3=nodes(n3,2);

    A = (0.5)*det([1 x1 y1; ...
                   1 x2 y2; ...
                   1 x3 y3]);

    m11 = (x2*y3 - x3*y2)/(2*A);
    m21 = (x3*y1 - x1*y3)/(2*A);
    m31 = (x1*y2 - y1*x2)/(2*A);
    m12 = (y2 - y3)/(2*A);
    m22 = (y3 - y1)/(2*A);
    m32 = (y1 - y2)/(2*A);
    m13 = (x3 - x2)/(2*A);
    m23 = (x1 - x3)/(2*A);
    m33 = (x2 -x1)/(2*A);

    B = [ m12  0  m22  0  m32  0; ...
           0  m13  0  m23  0  m33; ...
          m13 m12 m23 m22 m33 m32];

    k = thickness * A * B' * D * B ;

    K(2*n1-1,2*n1-1) = K(2*n1-1,2*n1-1) + k(1,1);
    K(2*n1-1,2*n1) = K(2*n1-1,2*n1) + k(1,2);
    K(2*n1-1,2*n2-1) = K(2*n1-1,2*n2-1) + k(1,3);
    K(2*n1-1,2*n2) = K(2*n1-1,2*n2) + k(1,4);
    K(2*n1-1,2*n3-1) = K(2*n1-1,2*n3-1) + k(1,5);
    K(2*n1-1,2*n3) = K(2*n1-1,2*n3) + k(1,6);
    K(2*n1,2*n1-1) = K(2*n1,2*n1-1) + k(2,1);
    K(2*n1,2*n1) = K(2*n1,2*n1) + k(2,2);
    K(2*n1,2*n2-1) = K(2*n1,2*n2-1) + k(2,3);
    K(2*n1,2*n2) = K(2*n1,2*n2) + k(2,4);
    K(2*n1,2*n3-1) = K(2*n1,2*n3-1) + k(2,5);
    K(2*n1,2*n3) = K(2*n1,2*n3) + k(2,6);
    K(2*n2-1,2*n1-1) = K(2*n2-1,2*n1-1) + k(3,1);
    K(2*n2-1,2*n1) = K(2*n2-1,2*n1) + k(3,2);
    K(2*n2-1,2*n2-1) = K(2*n2-1,2*n2-1) + k(3,3);
    K(2*n2-1,2*n2) = K(2*n2-1,2*n2) + k(3,4);
    K(2*n2-1,2*n3-1) = K(2*n2-1,2*n3-1) + k(3,5);
    K(2*n2-1,2*n3) = K(2*n2-1,2*n3) + k(3,6);
    K(2*n2,2*n1-1) = K(2*n2,2*n1-1) + k(4,1);
    K(2*n2,2*n1) = K(2*n2,2*n1) + k(4,2);
    K(2*n2,2*n2-1) = K(2*n2,2*n2-1) + k(4,3);
    K(2*n2,2*n2) = K(2*n2,2*n2) + k(4,4);
    K(2*n2,2*n3-1) = K(2*n2,2*n3-1) + k(4,5);
    K(2*n2,2*n3) = K(2*n2,2*n3) + k(4,6);
    K(2*n3-1,2*n1-1) = K(2*n3-1,2*n1-1) + k(5,1);
    K(2*n3-1,2*n1) = K(2*n3-1,2*n1) + k(5,2);
    K(2*n3-1,2*n2-1) = K(2*n3-1,2*n2-1) + k(5,3);
    K(2*n3-1,2*n2) = K(2*n3-1,2*n2) + k(5,4);
    K(2*n3-1,2*n3-1) = K(2*n3-1,2*n3-1) + k(5,5);
    K(2*n3-1,2*n3) = K(2*n3-1,2*n3) + k(5,6);
    K(2*n3,2*n1-1) = K(2*n3,2*n1-1) + k(6,1);
    K(2*n3,2*n1) = K(2*n3,2*n1) + k(6,2);
    K(2*n3,2*n2-1) = K(2*n3,2*n2-1) + k(6,3);
    K(2*n3,2*n2) = K(2*n3,2*n2) + k(6,4);
    K(2*n3,2*n3-1) = K(2*n3,2*n3-1) + k(6,5);
    K(2*n3,2*n3) = K(2*n3,2*n3) + k(6,6);
end

%% BOUNDARY CONDITIONS , DISPLACEMENTS %%

U = zeros(2*nnd,1);

cnt=0;

K_F = K;
F_F = F;

for i=1:size(bc,1)

    node = bc(i, 1);
    direction = bc(i, 2);
    value = bc(i, 3);

    U(2 * node - (2-direction)) = value;

    F_F = F_F - K_F(:,2*node-(2-direction)) * value;

    cnt=cnt+1;
    uu_zero(cnt)=2 * node - (2-direction);
end

K_F(:,uu_zero)=[];
K_F(uu_zero,:)=[];
F_F(uu_zero,:)=[];

uf = K_F \ F_F;

counter=1;

for i=1:2*nnd
    if i~=uu_zero
        U(i,1)=uf(counter,1);
        counter=counter+1;
    end
end


%% REACTION FORCES %%

reaction_forces = zeros(2*nnd,1);

for i=1:size(bc,1)

    node0 = bc(i, 1);
    direction0 = bc(i, 2);

    reaction_forces(2 * node0 - (2-direction0)) = K(2 * node0 - (2-direction0),:) * U;
end

%% STRESSES %%

SIGMA = zeros(3,nel);
EPS = zeros(3,nel);

for i=1:size(elements,1)

    n1=elements(i,1);
    n2=elements(i,2);
    n3=elements(i,3);

    x1=nodes(n1,1);
    x2=nodes(n2,1);
    x3=nodes(n3,1);

    y1=nodes(n1,2);
    y2=nodes(n2,2);
    y3=nodes(n3,2);

    A = (0.5)*det([1 x1 y1; ...
                   1 x2 y2; ...
                   1 x3 y3]);

    m11 = (x2*y3 - x3*y2)/(2*A);
    m21 = (x3*y1 - x1*y3)/(2*A);
    m31 = (x1*y2 - y1*x2)/(2*A);
    m12 = (y2 - y3)/(2*A);
    m22 = (y3 - y1)/(2*A);
    m32 = (y1 - y2)/(2*A);
    m13 = (x3 - x2)/(2*A);
    m23 = (x1 - x3)/(2*A);
    m33 = (x2 -x1)/(2*A);

    B = [ m12  0  m22  0  m32  0; ...
           0  m13  0  m23  0  m33; ...
          m13 m12 m23 m22 m33 m32];

    eps = B * [U(2*n1-1);
                 U(n1);
               U(2*n2-1);
                 U(n2);
               U(2*n3-1);
                 U(n3);];
    EPS(:,i) = eps;
    sigma = D * eps;
    SIGMA(:,i) = sigma;
end

%% DISPLAY RESULTS %%

disp('Nodal Displacements:');
for i = 1:nnd
    disp(['Node ' num2str(i) ' : X = ' num2str(U(2*i-1)) ', Y = ' num2str(U(2*i))]);
end

disp('Reaction Forces:');
for i = 1:nnd
    disp(['Node ' num2str(i) ' : X = ' num2str(reaction_forces(2*i-1)) ', Y = ' num2str(reaction_forces(2*i))]);
end

disp('Stresses:');
disp(SIGMA);

%% PLOTTING STRESSES CONTOURS %%

figure;
hold on;
title('Stresses Contours');

scatter(nodes(:,1), nodes(:,2), 'filled', 'MarkerFaceColor', 'b');
text(nodes(:,1), nodes(:,2), arrayfun(@num2str, 1:nnd, 'UniformOutput', false), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

for i = 1:size(elements,1)
    n1 = elements(i,1);
    n2 = elements(i,2);
    n3 = elements(i,3);
    plot([nodes(n1,1)+U(2*n1-1), nodes(n2,1)+U(2*n2-1), nodes(n3,1)+U(2*n3-1), nodes(n1,1)+U(2*n1-1)], ...
         [nodes(n1,2)+U(2*n1), nodes(n2,2)+U(2*n2), nodes(n3,2)+U(2*n3), nodes(n1,2)+U(2*n1)], 'k-');
end

for i = 1:size(elements,1)
    n1 = elements(i,1);
    n2 = elements(i,2);
    n3 = elements(i,3);
    
    s1 = SIGMA(1,i);
    s2 = SIGMA(2,i);
    s3 = SIGMA(3,i);

    mises = sqrt(s1^2 - s1*s2 + s2^2 + 3*s3^2);
    
    patch([nodes(n1,1)+U(2*n1-1), nodes(n2,1)+U(2*n2-1), nodes(n3,1)+U(2*n3-1), nodes(n1,1)+U(2*n1-1)], ...
          [nodes(n1,2)+U(2*n1), nodes(n2,2)+U(2*n2), nodes(n3,2)+U(2*n3), nodes(n1,2)+U(2*n1)], mises(1), ...
          'EdgeColor', 'k', 'LineWidth', 1.5);
end

colorbar;
xlabel('X-axis');
ylabel('Y-axis');
grid on;
