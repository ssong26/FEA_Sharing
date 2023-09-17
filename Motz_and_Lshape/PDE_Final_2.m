% =====================29-Apr-2018 Final Project============================
% Siyuan_Song Final Project for APMA 2560
% 29-Apr-Problem 1
%
% Program Description------------------------------------------------------
%
% 1. This program is to solve the uxx+uyy = -1 in L shape structure
%
% 2. The following variables will be calculated in the current program
%
% ==========================================================================
% =========================The Main Function================================
%
function PDE_Final_2()
clear;clc;close all
disp("FEM Analylsis Begin");
[xy_min,xy_max,xy_number] = deal(10,28,10);
xy_vector = floor(linspace(xy_min,xy_max,xy_number));
freedom_vector = zeros(1,xy_number);
L2_vector = zeros(1,xy_number);
Linf_vector = zeros(1,xy_number);
H1_vector = zeros(1,xy_number);
rh_vector = zeros(1,xy_number);
for i = 1 : xy_number
    %
    % ----------------------Part 1 Basic parameters----------------------------
    %
    % ---Structural Parameters
    [coeff,x_min,x_max,y_min,y_max,u_D] = deal(0,-1,1,-1,1,0);
    %
    %----------------------------Part 2 Generating Grid------------------------
    %
    % ----Basic parameters
    r_number = xy_vector(i);
    [x_number,y_number,alpha,layer] = deal(r_number,r_number,1.1,floor(1*r_number/2));
    % ----Uniform grid
    % [coord,connect,connect_number,nnode,rh,coord_pointer] = grid_generation(x_min,x_max,y_min,y_max,x_number,y_number);
    % ----Refine grid
    [coord,connect,connect_number,nnode,rh] = grid_iteractive(x_min,x_max,y_min,y_max,r_number,alpha,layer);
    %
    % ---------------------Part 3 Global Stiffness Matrix-----------------------
    %
    pressure = fem_Lshape(coord,connect,connect_number,nnode,u_D,coeff,rh);
    %
    % ========================Compare the result================================
    %
    % ================================
    % the contour of FEM
    % fem_contour(x_min,x_max,y_min,y_max,coord,pressure,connect)
    % ================================
    % the contour of exact solution
    % exact_contour(x_min,x_max,y_min,y_max,coord,nnode,connect)
    % ===============================
    % comparison figure of FEM and theory at the line y = y_line
    % ===============================
    % stress_Factor
    % stress_factor2(nnode,coord,pressure,rh)
    % ===============================
    % write output
    % fem_output(coord,pressure,connect)
    % ===============================
    % error contour
    % err_contour(coord,connect,pressure,nnode)
    % derr_contour(coord,connect,connect_number,pressure)
    % ===============================
    % error analysis
    [L2_vector(i),Linf_vector(i),H1_vector(i)] = error_norm(coord,connect,connect_number,pressure);
    freedom_vector(i) = nnode;
    rh_vector(i) = rh;
    %===================================end====================================
end
figure;
plot(log10(freedom_vector),log10(L2_vector));hold on
plot(log10(freedom_vector),log10(Linf_vector));hold on
plot(log10(freedom_vector),log10(H1_vector));
legend('L2 norm','Linf norm','H1 norm')
xlabel('log10(freedom)');ylabel('log10(error)')
% figure
% plot(log10(freedom_vector(2:xy_number)),-2 * (log10(L2_vector(2:xy_number))-log10(L2_vector(1:(xy_number-1)))) ./ (log10(freedom_vector(2:xy_number))-log10(freedom_vector(1:(xy_number-1)))));hold on
% plot(log10(freedom_vector(2:xy_number)),-2 * (log10(Linf_vector(2:xy_number))-log10(Linf_vector(1:(xy_number-1)))) ./ (log10(freedom_vector(2:xy_number))-log10(freedom_vector(1:(xy_number-1)))));hold on
% plot(log10(freedom_vector(2:xy_number)),-2 * (log10(H1_vector(2:xy_number))-log10(H1_vector(1:(xy_number-1)))) ./ (log10(freedom_vector(2:xy_number))-log10(freedom_vector(1:(xy_number-1)))));hold on
% legend('L2 norm','Linf norm','H1 norm')
% xlabel('log10(freedom)');ylabel('Error order')
figure;
plot(log10(freedom_vector(2:xy_number)),(log10(L2_vector(2:xy_number))-log10(L2_vector(1:(xy_number-1))))./(log10(rh_vector(2:xy_number))-log10(rh_vector(1:(xy_number-1)))));hold on
plot(log10(freedom_vector(2:xy_number)),(log10(Linf_vector(2:xy_number))-log10(Linf_vector(1:(xy_number-1))))./(log10(rh_vector(2:xy_number))-log10(rh_vector(1:(xy_number-1)))));hold on
plot(log10(freedom_vector(2:xy_number)),(log10(H1_vector(2:xy_number))-log10(H1_vector(1:(xy_number-1))))./(log10(rh_vector(2:xy_number))-log10(rh_vector(1:(xy_number-1)))));hold on
legend('L2 norm','Linf norm','H1 norm')
xlabel('log10(freedom)');ylabel('Error order')
end
%=================Some Functions may be used in the Main function==========
% ========================= Grid Generation
% Method 1 : generate the uniform grid (wrong)
function [coord,connect,connect_number,nnode,rh,coord_pointer] = grid_generation(x_min,x_max,y_min,y_max,x_number,y_number)
disp("Generate Uniform Grid")
nnode = (2 * x_number - 1) * y_number + x_number * (y_number - 1);    % number of node in the structure
coord = zeros(nnode,2);        % coordinates of the node
coord_pointer = zeros((2 * y_number - 1),(2 * x_number - 1));
% ---Generate the node and the basis connect
% from bottom to upper part
x_temp_bottom = linspace(x_min,x_max,(2 * x_number - 1));     % record the x position for bottom part
y_temp_bottom = linspace(y_min,0,y_number);     % record the y position for bottom part
rh = (0-y_min)/(y_number-1);
% bottom record
for i = 1 : y_number
    for j = 1 : (2 * x_number - 1)
        node_count = (i - 1) * (2 * x_number - 1) + j;
        coord(node_count,1)= x_temp_bottom(j);
        coord(node_count,2)= y_temp_bottom(i);
        coord_pointer(i,j) = node_count;
    end
end
x_temp_upper = linspace(x_min,0,x_number);     % record the x position for upper part
y_temp_upper = linspace((0-y_min)/(y_number-1),y_max,y_number-1);     % record the y position for upper part
% upper record
for i = 1 : (y_number-1)
    for j = 1 : x_number
        node_count = y_number * (2 * x_number - 1)+ (i - 1) * x_number + j;
        coord(node_count,1)= x_temp_upper(j);
        coord(node_count,2)= y_temp_upper(i);
        coord_pointer(i + y_number,j) = node_count;
    end
end
connect_temp = delaunay(coord(:,1),coord(:,2));%generate the connectivity array
% delete the overlap triangle
connect = [];
connect_temp_count = 1;
for i = 1 : size(connect_temp(:,1))
    % midpoints
    triangle_x = coord(connect_temp(i,1),1) + coord(connect_temp(i,2),1) + coord(connect_temp(i,3),1);
    triangle_y = coord(connect_temp(i,1),2) + coord(connect_temp(i,2),2) + coord(connect_temp(i,3),2);
    if triangle_x > 0 && triangle_y > 0
    else
        connect(connect_temp_count,1) = connect_temp(i,1);
        connect(connect_temp_count,2) = connect_temp(i,2);
        connect(connect_temp_count,3) = connect_temp(i,3);
        connect_temp_count = connect_temp_count + 1;
    end
end
connect_number=size(connect,1);
% figure;
% triplot(connect,coord(:,1),coord(:,2),'r');
end
% Method 2 : generate the refining grid (wrong)
function [coord,connect,connect_number,nnode,rh] = grid_iteractive(x_min,x_max,y_min,y_max,r_number,alpha,layer)
disp("Generate Local Refining Grid");
coord = [];
number = 2;
coord(1,1) = 0; coord(1,2) = 0;
rh = 1/((alpha^layer-1)/(alpha-1)+(r_number-1-layer)*alpha^(layer-1));
% obtain r_vector
r_vector = zeros(r_number,1);
for i = 1 : layer
    r_vector(i+1) = r_vector(i) + rh * alpha^(i-1);
end
for i = (layer + 2) : r_number
    r_vector(i) = r_vector(i-1) + rh * alpha^(layer-1);
end
r_distance = zeros(r_number-1,1);
for i = 1 : layer
    r_distance(i) = rh * alpha^(i-1);
end
for i = (layer + 1) : (r_number-1)
    r_distance(i) = rh * alpha^(layer-1);
end
r_vector(r_number) = 1;
% loop over all point in r-direction
for i = 1 : (r_number-1)
    increment = r_distance(i);
    temp_number = ceil(r_vector(i+1)/increment);
    increment = r_vector(i+1)/temp_number;
    %
    coord(number,1) = 0 ;
    coord(number,2) = r_vector(i+1);
    number = number + 1;
    r_temp = linspace(increment,r_vector(i+1),temp_number);
    for j = 1 : temp_number
        coord(number,1) = - r_temp(j);
        coord(number,2) = r_vector(i+1);
        number = number + 1;
    end
    for j = 1 : temp_number
        coord(number,1) = - r_vector(i+1);
        coord(number,2) = r_vector(i+1) - r_temp(j);
        number = number + 1;
    end
    for j = 1 : temp_number
        coord(number,1) = - r_vector(i+1);
        coord(number,2) = - r_temp(j);
        number = number + 1;
    end
    for j = 1 : temp_number
        coord(number,1) = - r_vector(i+1) + r_temp(j);
        coord(number,2) = - r_vector(i+1);
        number = number + 1;
    end
    for j = 1 : temp_number
        coord(number,1) = 0 + r_temp(j);
        coord(number,2) = - r_vector(i+1);
        number = number + 1;
    end
    for j = 1 : temp_number
        coord(number,1) =   r_vector(i+1);
        coord(number,2) = - r_vector(i+1) + r_temp(j);
        number = number + 1;
    end
end
%
nnode = length(coord(:,1));
connect_temp = delaunayTriangulation(coord(:,1),coord(:,2));%generate the connectivity array
% delete the overlap triangle
connect = [];
connect_temp_count = 1;
for i = 1 : size(connect_temp(:,1))
    % midpoints
    triangle_x = coord(connect_temp(i,1),1) + coord(connect_temp(i,2),1) + coord(connect_temp(i,3),1);
    triangle_y = coord(connect_temp(i,1),2) + coord(connect_temp(i,2),2) + coord(connect_temp(i,3),2);
    if triangle_x > 0 && triangle_y > 0
    else
        connect(connect_temp_count,1) = connect_temp(i,1);
        connect(connect_temp_count,2) = connect_temp(i,2);
        connect(connect_temp_count,3) = connect_temp(i,3);
        connect_temp_count = connect_temp_count + 1;
    end
end
connect_number=size(connect,1);
% figure;
% triplot(connect,coord(:,1),coord(:,2),'r');
rh = rh * alpha^(layer-1);
end
% Find (x,y) belong to which element (just fit this problem uniform grid)
function node_pointer = element_find(x,y,coord,coord_pointer,r_number)
rh = 1 / (r_number - 1);
x0 = floor(x/rh);
y0 = floor(y/rh);
node_pointer = zeros(3,1);
judgement = false;
% there are a lot of discussion for different occassions
% Case one: x = y = 0
if x == 0 && y == 0 && judgement == false
    node_pointer(1) = coord_pointer(r_number  ,r_number  );
    node_pointer(2) = coord_pointer(r_number-1,r_number  );
    node_pointer(3) = coord_pointer(r_number  ,r_number-1);
    judgement = true;
end
% Case two: x = 0, 0 < y < 1; 
if x == 0 && y * (1 - y) > 0 && judgement == false
    node_pointer(1) = coord_pointer(r_number + y0,     r_number);
    node_pointer(2) = coord_pointer(r_number + y0 + 1, r_number);
    node_pointer(3) = coord_pointer(r_number + y0 + 1, r_number-1);
    judgement = true;
end
% Case three: x = 0, y = 1; 
if x == 0 && y == 1 && judgement == false
    node_pointer(1) = coord_pointer(2 * r_number - 1, r_number);
    node_pointer(2) = coord_pointer(2 * r_number - 2, r_number);
    node_pointer(3) = coord_pointer(2 * r_number - 1, r_number - 1);
    judgement = true;
end
% Case four: -1 < x < 0, y = 1; 
if x > -1 && x < 0 && y == 1 && judgement == false
    node_pointer(1) = coord_pointer(2 * r_number - 1, r_number + x0);
    node_pointer(2) = coord_pointer(2 * r_number - 1, r_number + x0 + 1);
    node_pointer(3) = coord_pointer(2 * r_number - 2, r_number + x0 + 1);
    judgement = true;
end
% Case five: x = -1, y = 1; 
if x == -1 && y == 1 && judgement == false
    node_pointer(1) = coord_pointer(2 * r_number - 1, 1);
    node_pointer(2) = coord_pointer(2 * r_number - 1, 2);
    node_pointer(3) = coord_pointer(2 * r_number - 2, 2);
    judgement = true;
end
% Case six: x = -1, -1 < y < 1;
if x == -1 && y > -1 && y < 1 && judgement == false
    node_pointer(1) = coord_pointer(r_number + y0, 1);
    node_pointer(2) = coord_pointer(r_number + y0 + 1, 1);
    node_pointer(3) = coord_pointer(r_number + y0, 2);
    judgement = true;
end
% Case seven: x = -1, y = -1;
if x == -1 && y == -1 && judgement == false
    node_pointer(1) = coord_pointer(1, 1);
    node_pointer(2) = coord_pointer(2, 1);
    node_pointer(3) = coord_pointer(1, 2);
    judgement = true;
end
% Case eight: 1 > x > -1, y = -1;
if x > -1 && x < 1 && y == -1 && judgement == false
    node_pointer(1) = coord_pointer(1, r_number + x0);
    node_pointer(2) = coord_pointer(2, r_number + x0);
    node_pointer(3) = coord_pointer(1, r_number + x0 + 1);
    judgement = true;
end
% Case nine: x = 1, y = -1;
if x == 1 && y == -1 && judgement == false
    node_pointer(1) = coord_pointer(1, 2 * r_number - 1);
    node_pointer(2) = coord_pointer(2, 2 * r_number - 1);
    node_pointer(3) = coord_pointer(2, 2 * r_number - 2);
    judgement = true;
end
% Case ten: x = 1, 0 > y > -1;
if x == 1 && y > -1 && y < 0 && judgement == false
    node_pointer(1) = coord_pointer(r_number + y0, 2 * r_number - 1);
    node_pointer(2) = coord_pointer(r_number + y0 + 1, 2 * r_number - 1);
    node_pointer(3) = coord_pointer(r_number + y0 + 1, 2 * r_number - 2);
    judgement = true;
end
% Case eleven: x = 1, y = 0;
if x == 1 && y == 0 && judgement == false
    node_pointer(1) = coord_pointer(r_number, 2 * r_number - 1);
    node_pointer(2) = coord_pointer(r_number, 2 * r_number - 2);
    node_pointer(3) = coord_pointer(r_number - 1, 2 * r_number - 1);
    judgement = true;
end
% Case twelve: 0 < x < 1, y = 0;
if x < 1 && x > 0 && y == 0 && judgement == false
    node_pointer(1) = coord_pointer(r_number, r_number + x0);
    node_pointer(2) = coord_pointer(r_number, r_number + x0 + 1);
    node_pointer(3) = coord_pointer(r_number - 1, r_number + x0 + 1);
    judgement = true;
end
% Case last: else
if judgement == false
    if ((x - x0*rh) + (y - y0*rh)) <= rh
        node_pointer(1) = coord_pointer(r_number + y0, r_number + x0);
        node_pointer(2) = coord_pointer(r_number + y0 + 1, r_number + x0);
        node_pointer(3) = coord_pointer(r_number + y0, r_number + x0 + 1);
        judgement = true;
    else
        node_pointer(1) = coord_pointer(r_number + y0 + 1, r_number + x0 + 1);
        node_pointer(2) = coord_pointer(r_number + y0 + 1, r_number + x0);
        node_pointer(3) = coord_pointer(r_number + y0, r_number + x0 + 1);
        judgement = true;
    end
end
% for test
% figure;
% scatter(x,y,'*');hold on;
% scatter(coord(node_pointer(1),1),coord(node_pointer(1),2));hold on;
% scatter(coord(node_pointer(2),1),coord(node_pointer(2),2));hold on;
% scatter(coord(node_pointer(3),1),coord(node_pointer(3),2));hold on;
end
% ========================= FEM Calculation
% finite element program
function pressure = fem_Lshape(coord,connect,connect_number,nnode,u_D,coeff,rh)
disp("Stiffness Matrix Assemble...");
Stif=zeros(nnode,nnode);
for lmn=1:connect_number    % Loop over all the elements
    %
    %   Set up the stiffness for the current element
    %
    a = connect(lmn,1);
    b = connect(lmn,2);
    c = connect(lmn,3);
    coord_element = zeros(3,2);
    coord_element(1,:)=coord(a,:);
    coord_element(2,:)=coord(b,:);
    coord_element(3,:)=coord(c,:);
    n_nodes=3;n_points=3;
    kel = elstif(coord_element,coeff,n_nodes,n_points);
    %
    %   Add the current element stiffness to the global stiffness
    %
    for i = 1 : 3
        for j = 1 : 3
            Stif(connect(lmn,i),connect(lmn,j)) = Stif(connect(lmn,i),connect(lmn,j)) + kel(i,j);
        end
    end
end
% figure
% spy(Stif)
% p = symrcm(Stif); 
% figure
% spy(Stif(p,p))
%
%================
%
% ==========================================================================
% =========================Boundary Condition 1=============================
%
% fixnode(m,1) store the Serial Number of the node on the dirichlet boundary 
% fixnode(m,2) store the value of the node on the dirichlet boundary.
% There is a trouble that, whether or not I need to include the point (0,0) into the boundary condition.
% At present, I determine to include (0,0) in the boundary condition.
disp("Boundary Conditions Setting...");
fixnodes_number = 1;
fixnodes = [];
fixnodes_count = 1;
for i = 1 : nnode
    x_temp = coord(i,1);
    y_temp = coord(i,2);
    judgement = 0;
    if (x_temp - 1)*(x_temp + 1)*(y_temp + 1)*(y_temp - 1) <= 10^(-10) && judgement == 0
        fixnodes(fixnodes_count,1) = i;
        fixnodes(fixnodes_count,2) = u_D;
        fixnodes_count = fixnodes_count + 1;
        judgement = 1;
    end
    if abs(x_temp) <= min(10^(-20),abs(rh/10)) && y_temp*(1-y_temp) >= min(10^(-20),abs(rh/10)) && judgement == 0
        fixnodes(fixnodes_count,1) = i;
        fixnodes(fixnodes_count,2) = u_D;
        fixnodes_count = fixnodes_count + 1;
        judgement = 1;
    end
    if abs(y_temp) <= min(10^(-20),abs(rh/10)) && x_temp*(1-x_temp) >= min(10^(-20),abs(rh/10)) && judgement == 0
        fixnodes(fixnodes_count,1) = i;
        fixnodes(fixnodes_count,2) = u_D;
        fixnodes_count = fixnodes_count + 1;
        judgement = 1;
    end
    if abs(x_temp) <= min(10^(-20),abs(rh/10))  && abs(y_temp) <= min(10^(-20),abs(rh/10)) && judgement == 0
        fixnodes(fixnodes_count,1) = i;
        fixnodes(fixnodes_count,2) = u_D;
        fixnodes_count = fixnodes_count + 1;
    end
end
fixnodes_number = fixnodes_count - 1;
% figure
% scatter(coord(fixnodes(:,1),1),coord(fixnodes(:,1),2));
%
% ==================== Assemble the r vector ==============================
%
disp("Residual Setting...");
resid = zeros(nnode,1);
for lmn = 1 : connect_number    % loop over all elements
    %
    %   Set up the residual for the current element
    %
    a = connect(lmn,1);
    b = connect(lmn,2);
    c = connect(lmn,3);
    coord_element = zeros(3,2);
    coord_element(1,:)=coord(a,:);
    coord_element(2,:)=coord(b,:);
    coord_element(3,:)=coord(c,:);
    n_nodes=3;n_points=3;
    rel = residual(coord_element,n_nodes,n_points);
    %
    %   Add the current element stiffness to the global stiffness
    %
    for i = 1 : 3
        resid(connect(lmn,i)) = resid(connect(lmn,i)) + rel(i);
    end
end
%
% ==================== Assemble Boundary condition ========================
%
% Modify the global stiffness and residual to include constraints
% The symmetry must be satisfied
%
for i=1:fixnodes_number
    for j=1:nnode
        resid(j) = resid(j) - Stif(j,fixnodes(i,1)) * fixnodes(i,2);
        Stif(fixnodes(i,1),j) = 0;
        Stif(j,fixnodes(i,1)) = 0;
    end
    resid(fixnodes(i,1)) = fixnodes(i,2);
    Stif(fixnodes(i,1),fixnodes(i,1)) = 1;
end
% ================== Solve the FEM equations ==============================
%
disp("Getting The Result...");
pressure = Stif\resid;
end
% ========================= Result Discussion
% contour of the fem solutions
function fem_contour(x_min,x_max,y_min,y_max,coord,pressure,connect)
%{
x_mesh_number = 100;
y_mesh_number = 100;
xq = linspace(x_min,x_max,x_mesh_number);yq = linspace(y_min,y_max,y_mesh_number);
[Xq,Yq] = meshgrid(xq,yq);
Zq=griddata(coord(:,1),coord(:,2),pressure,Xq,Yq,'v4');
figure;
pcolor(Xq,Yq,Zq);
colorbar;
%}
disp("generate numerical contour...");
figure
trimesh(connect,coord(:,1),coord(:,2),pressure); % hand in this plot
xlabel('x'), ylabel('y'), zlabel('u(x,y)'), az = 70; el = 60; view(az, el);
end
% contour of the derivatives
function derr_contour(coord,connect,connect_number,pressure)
ux_error = zeros(connect_number,2);
ux_coord = zeros(connect_number,2);
for lmn = 1 : connect_number
    % get the position of the three nodes of the triangle element
    x1_node = coord(connect(lmn,1:3),1);
    x2_node = coord(connect(lmn,1:3),2);
    u_node = pressure(connect(lmn,1:3));
    % obtain the position of the center points
    x1_center = sum(x1_node)/3;
    x2_center = sum(x2_node)/3;
    % numerical results
    [~,u_x1,u_x2] = u_predict(x1_center,x2_center,x1_node,x2_node,u_node);
    % theoretical results
    %{
    dx = 10^(-5);
    u_center_Motz = L_shape(x1_center,x2_center);
    u_center_Motz_x1 = (L_shape(x1_center - dx, x2_center + 0 ) - u_center_Motz)/(-dx);
    u_center_Motz_x2 = (L_shape(x1_center + 0 , x2_center - dx) - u_center_Motz)/(-dx);
    %}
    u_temp = L_exact([x1_center,x2_center]);
    u_center_Motz = u_temp(1);
    u_center_Motz_x1 = u_temp(2);
    u_center_Motz_x2 = u_temp(3);
    % 
    ux_error(lmn,1) = abs(u_x1- u_center_Motz_x1);
    ux_error(lmn,2) = abs(u_x2- u_center_Motz_x2);
    % ux_error(lmn) = sqrt( u_center_Motz_x1^2 + u_center_Motz_x2^2 );
    ux_coord(lmn,1) = x1_center;
    ux_coord(lmn,2) = x2_center;
end
figure
trimesh(connect,ux_coord(:,1),ux_coord(:,2),ux_error(:,1)); % hand in this plot
xlabel('x'), ylabel('y'), zlabel('u(x,y)'), az = 70; el = 60; view(az, el);
figure
trimesh(connect,ux_coord(:,1),ux_coord(:,2),ux_error(:,2)); % hand in this plot
xlabel('x'), ylabel('y'), zlabel('u(x,y)'), az = 70; el = 60; view(az, el);
end
% contour of the exact solutions
function exact_contour(x_min,x_max,y_min,y_max,coord,nnode,connect)
%{
x_mesh_number = 100;
y_mesh_number = 100;
xq = linspace(x_min,x_max,x_mesh_number);yq = linspace(y_min,y_max,y_mesh_number);
[Xq,Yq] = meshgrid(xq,yq);
pressure_exact = zeros(nnode,1);
for i = 1 : nnode
    pressure_exact(i) = L_shape(coord(i,1),coord(i,2));
end
Zq=griddata(coord(:,1),coord(:,2),pressure_exact,Xq,Yq,'v4');
figure;
gca = pcolor(Xq,Yq,Zq);
colorbar;
set(gca,'LineStyle','none');
%}
pressure_exact = zeros(nnode,1);
for i = 1 : nnode
%    pressure_exact(i) = L_shape(coord(i,1),coord(i,2));
    u_temp = L_exact(coord(i,:));
    pressure_exact(i) = u_temp(1);
end
disp("generate exact contour...");
figure
trimesh(connect,coord(:,1),coord(:,2),pressure_exact); % hand in this plot
xlabel('x'), ylabel('y'), zlabel('u(x,y)'), az = 70; el = 60; view(az, el);
end
% compare the value of u on the position y = y_line
function compare_yline(x_min,x_max,y_min,y_max,x_number,y_number,coord,pressure,y_yline)
figure;
% theory
y_yline = -0.2;
x_yline_number = 200;
x_yline = linspace(x_min,x_max,x_yline_number);
u_yline_exact = zeros(1,x_yline_number);
for i = 1 : x_yline_number
    u_yline_exact(i) = L_shape(x_yline(i),y_yline);
end
plot(x_yline,u_yline_exact,'LineWidth',4);hold on
% fem
x_temp_bottom = linspace(x_min,x_max,(2 * x_number - 1));
y_position = floor((y_yline + y_max)/abs(y_min)*y_number) + 1; 
scatter(x_temp_bottom,pressure(((y_position - 1)*(2*x_number - 1)+1):y_position*(2*x_number - 1)),150)
legend('theoretical prediction','numerical simulation')
xlabel('X-direction');ylabel('u(x,y)')
end
% Method 1 : Calculate the stress concentration factor via the min point
function stress_factor(nnode,coord,pressure,rh)
point_r = 1000;
A1 = 1;
A1_error = abs(0.40192487 - A1)/0.40192487;
for i = 1:nnode
    x = coord(i,1); y = coord(i,2);
    r = sqrt(x^2 + y^2);
    if r < point_r && r > 0
        if abs(x)>10^(-10) && abs(y)>10^(-10)
            if x < 0 && y > 0
                theta = pi - atan(-y/x);
            end
            if x < 0 && y < 0
                theta = pi + atan(y/x);
            end
            if x > 0 && y < 0
                theta = 2*pi + atan(y/x);
            end
            up = -r^2 / (6*pi) * (3*pi/2 - 2*log(r)*sin(2*theta) - (2*theta+3*pi/2)*cos(2*theta));
            A1 = ( pressure(i) - up ) / sin(2*theta/3-pi/3) / r^(2/3);
            A1_error = abs(0.40192487 - A1)/0.40192487;
            point_r = r;
        end
    end
end
disp('r = ')
disp(point_r);
disp('A1 = ')
disp(A1);
disp('A1_error = ')
disp(A1_error);
end
% Method 2 : Calculate the stress concentration factor via several points
function stress_factor2(nnode,coord,pressure,rh)
A1 = 1;
A1_error = abs(0.40192487 - A1)/0.40192487;
s = 0; number = 0;
for i = 1:nnode
    x = coord(i,1); y = coord(i,2);
    r = sqrt(x^2 + y^2);
    if r < abs(0.2) && r > 0
        if abs(x)>0 && abs(y)>0
            if x < 0 && y > 0
                theta = pi - atan(-y/x);
            end
            if x < 0 && y < 0
                theta = pi + atan(y/x);
            end
            if x > 0 && y < 0
                theta = 2*pi + atan(y/x);
            end
            up = -r^2 / (6*pi) * (3*pi/2 - 2*log(r)*sin(2*theta) - (2*theta+3*pi/2)*cos(2*theta));
            s = s + ( pressure(i) - up ) / sin(2*theta/3-pi/3) / r^(2/3);
            number = number + 1;
        end
    end
end
A1 = s / number;
A1_error = abs(0.40192487 - A1)/0.40192487;
disp('A1 = ')
disp(A1);
disp('A1_error = ')
disp(A1_error);
end
% ========================= Subfunction for FEM
% integrationpoints function for 2D, we can obtain the position and weight
function xi = abq_UEL_2D_integrationpoints(n_points, n_nodes)
xi=zeros(n_nodes,3);        % xi=(x,y,w)
if n_nodes==3               % nodes of the elements
    if n_points==3          % number of integration points
        xi(1,1)=0.5;xi(1,2)=0.5;xi(1,3)=1/6; % integration point 1
        xi(2,1)=0.0;xi(2,2)=0.5;xi(2,3)=1/6; % integration point 2
        xi(3,1)=0.5;xi(3,2)=0.0;xi(3,3)=1/6; % integration point 3
    end
end
end
% shapefunction for 2D, we can obtain the N and dNdx.
function f = abq_UEL_2D_shapefunctions(xi,n_points,n_nodes)
f=zeros(n_nodes,3);
if n_nodes==3
    if n_points==3
        f(1,1)=xi(1);f(1,2)=1;f(1,3)=0;
        f(2,1)=xi(2);f(2,2)=0;f(2,3)=1;
        f(3,1)=1-xi(1)-xi(2);f(3,2)=-1;f(3,3)=-1;
    end
end
end
% stiffness function, we can obtain the stiffness matrix==
function kel = elstif(coord,coeff,n_nodes,n_points)
xi = abq_UEL_2D_integrationpoints(n_points, n_nodes);
kel = 0;
for i = 1 : 3 % integration over all points
    f = abq_UEL_2D_shapefunctions(xi(i,1:2),n_points,n_nodes);
    %
    dNdxi = zeros(3,2);dNdxi(:,1) = (f(:,2))';dNdxi(:,2) = (f(:,3))';
    %
    dxdxi = coord'*dNdxi;
    dxidx = inv(dxdxi);
    %
    dNdx = dNdxi*dxidx;
    determinate = dxdxi(1,1)*dxdxi(2,2)-dxdxi(1,2)*dxdxi(2,1);
    % ---C matrix
    if n_points ==3
        if n_nodes==3
            B_wave=zeros(3,2);
            B_wave(:,1) = dNdx(:,1);B_wave(:,2) = dNdx(:,2);
        end
    end
    C = B_wave*B_wave'*xi(i,3)*determinate;
    % ---T matrixv
    if n_points == 3
        if n_nodes == 3
            B_line=zeros(3,1);
            B_line(:,1)=f(:,1);
        end
    end
    T = B_line*B_line'*xi(i,3)*determinate;
    kel = kel + 1*(C-coeff^2*T);
end
end
% calculate the residual of the fem
function rel = residual(coord,n_nodes,n_points)
xi = abq_UEL_2D_integrationpoints(n_points, n_nodes);
rel = 0;
for i = 1 : 3 % integration over all points
    f = abq_UEL_2D_shapefunctions(xi(i,1:2),n_points,n_nodes);
    %
    dNdxi = zeros(3,2);dNdxi(:,1) = (f(:,2))';dNdxi(:,2) = (f(:,3))';
    %
    dxdxi = coord'*dNdxi;
    dxidx = inv(dxdxi);
    %
    dNdx = dNdxi*dxidx;
    determinate = dxdxi(1,1)*dxdxi(2,2)-dxdxi(1,2)*dxdxi(2,1);
    % ---r vector
    if n_points == 3
        if n_nodes == 3
            B_line=zeros(3,1);
            B_line(:,1)=f(:,1);
        end
    end
    T = B_line*xi(i,3)*determinate;
    rel = rel + T;
end
end
% ========================= Exact Solution and Writing Output
% To fit Ak
function ak = ak_fit()
clear;clc
coord_stru = load('coord.mat');
pressure_stru = load('pressure.mat');
coord = coord_stru.coord;
pressure = pressure_stru.pressure;
nnode = length(coord(:,1));
fit_number = 20;
y_data = zeros(nnode,1);
x_data = zeros(nnode,fit_number);
for i = 1:nnode
    x = coord(i,1); y = coord(i,2);
    r = sqrt(x^2+y^2);
    % First coordinate
    if x > 0 && y == 0
        theta = 2*pi;
    end
    if x == 0 && y > 0
        theta = pi/2;
    end
    % Second coordinate
    if x < 0 && y > 0
        theta = pi - atan(-y/x);
    end
    if x < 0 && y == 0
        theta = pi;
    end
    % Third coordinate
    if x < 0 && y < 0
        theta = pi + atan(y/x);
    end
    if x == 0 && y < 0
        theta = 3/2 * pi;
    end
    % Fourth coordinate
    if x > 0 && y < 0
        theta = 2*pi + atan(y/x);
    end
    % second, obtain the real solution
    if r > 0
        up = -r^2 / (6*pi) * (3*pi/2 - 2*log(r)*sin(2*theta) - (2*theta - 5*pi/2)*cos(2*theta));
        y_data(i) = pressure(i) - up;
        for j = 1:fit_number
            x_data(i,j) = r^((2*j-1) * 2/3) * sin((2*j-1)*(2*theta-pi)/3);
        end
    end
end
ak = regress(reshape(y_data,nnode,1),x_data);
end
% Exact solution function by given form (very useful)
function s = L_shape(x,y)
% initially, determine the coefficient
coeff = [      0.401924749945678
   0.093647322392899
  -0.009382266873999
  -0.029886646782791
  -0.008358994236376
  -0.004727929057437
  -0.001541324424179
  -0.001092134923554
  -0.000697767541724
  -0.000538726404143
  -0.000366251120773
  -0.000232513632145
  -0.000154344089805
  -0.000112407286977
  -0.000051702524179
  -0.000033989886018
  -0.000023728176110
  -0.000005727204713
  -0.000003045494283
  -0.000002103173041];
%first, obtain the radial coordinates
r = sqrt(x^2+y^2);
theta = 0;
% First coordinate
if x > 0 && y == 0
    theta = 2*pi;
end
if x > 0 && y > 0 && y > x
    theta = atan(y/x);
end
if x > 0 && y > 0 && x <= y
    theta = 2*pi + atan(y/x);
end
if x == 0 && y > 0
    theta = pi/2;
end 
% Second coordinate
if x < 0 && y > 0
    theta = pi - atan(-y/x);
end
if x < 0 && y == 0
    theta = pi;
end
% Third coordinate
if x < 0 && y < 0
    theta = pi + atan(y/x);
end
if x == 0 && y < 0
    theta = 3/2 * pi;
end
% Fourth coordinate
if x > 0 && y < 0
    theta = 2*pi + atan(y/x);
end
% second, obtain the real solution
if r > 0
    uh = 0;
    for j = 1:20
        uh = uh + coeff(j) * r^((2*j-1) * 2/3) * sin((2*j-1)*(2*theta-pi)/3);
    end
    up = -r^2 / (6*pi) * (3*pi/2 - 2*log(r)*sin(2*theta) - (2*theta - 5*pi/2)*cos(2*theta));
    s = uh + up;
end
if r == 0
    s = 0;
end
end
% Writing output for the current problem
function fem_output(coord,pressure,connect)
nnode = length(coord(:,1));
node_value = zeros(nnode,3);
for i = 1:nnode
    node_value(i,1) = coord(i,1);
    node_value(i,2) = coord(i,2);
    node_value(i,3) = pressure(i);
end
disp("Write Output...");
txt_write('node_value.txt',node_value);
txt_write('connect.txt',connect);
disp("Write Successfully...");
end
% Basic txt write program
function txt_write(name,matrix)
a = matrix;
[m,n] = size(a);
fid = fopen(name,'w');
for i=1:1:m
    for j=1:1:n
        if j==n
            fprintf(fid,'%g\n',a(i,j));
        else
            fprintf(fid,'%g\t',a(i,j));
        end
    end
end
fclose(fid);
end
% Obtain the exact solution. (very slow)
function u_exact = L_exact(coord_exact)
% read data, and save to [coord,connect,pressure]
% connect_stru = load('connect.mat');
coord_stru = load('coord.mat');
pressure_stru = load('pressure.mat');
coord_pointer_stru = load('coord_pointer.mat');
% connect = connect_stru.connect;
coord = coord_stru.coord;
pressure = pressure_stru.pressure;
coord_pointer = coord_pointer_stru.coord_pointer;
% connect_number = length(connect(:,1));
r_number = (length(coord_pointer(:,1)) + 1)/2;
% begin
coord_exact_number = length(coord_exact(:,1));
u_exact = zeros(coord_exact_number,3); % u, ux, uy
for j = 1 : coord_exact_number % loop over all integration points
    node_pointer = element_find(coord_exact(j,1),coord_exact(j,2),coord,coord_pointer,r_number);
    coord_element = zeros(3,2);
    pressure_element = zeros(3,1);
    for ttt = 1 : 3
        coord_element(ttt,1) = coord(node_pointer(ttt),1);
        coord_element(ttt,2) = coord(node_pointer(ttt),2);
        pressure_element(ttt) = pressure(node_pointer(ttt));
    end
    [u_exact(j,1),u_exact(j,2),u_exact(j,3)] = u_predict(coord_exact(j,1),coord_exact(j,2),coord_element(:,1),coord_element(:,2),pressure_element);
end
end
% To test if the point coord_p in the element coord_element
function judgement = element_judge(coord_p,coord_element)
judgement = false;
coord_a = coord_element(1,:);
coord_b = coord_element(2,:);
coord_c = coord_element(3,:);
v0 = reshape(coord_c,2,1) - reshape(coord_a,2,1);
v1 = reshape(coord_b,2,1) - reshape(coord_a,2,1);
v2 = reshape(coord_p,2,1) - reshape(coord_a,2,1);
Stif = [v0(1),v1(1);
     v0(2),v1(2)];
Resid = [v2(1);v2(2)];
uv = Stif\Resid;
if uv(1) >= 0 && uv(2) >= 0 && sum(uv) <= 1
    judgement = true;
end
end
% To obtain the predict solution at the point x1,x2
function [u,u_x1,u_x2] = u_predict(x1,x2,x1_node,x2_node,u_node)
Na = ( (x2-x2_node(2))*(x1_node(3)-x1_node(2)) - (x1-x1_node(2))*(x2_node(3)-x2_node(2)) ) / ( (x2_node(1)-x2_node(2))*(x1_node(3)-x1_node(2)) - (x1_node(1)-x1_node(2))*(x2_node(3)-x2_node(2)) );
Nb = ( (x2-x2_node(3))*(x1_node(1)-x1_node(3)) - (x1-x1_node(3))*(x2_node(1)-x2_node(3)) ) / ( (x2_node(2)-x2_node(3))*(x1_node(1)-x1_node(3)) - (x1_node(2)-x1_node(3))*(x2_node(1)-x2_node(3)) );
Nc = ( (x2-x2_node(1))*(x1_node(2)-x1_node(1)) - (x1-x1_node(1))*(x2_node(2)-x2_node(1)) ) / ( (x2_node(3)-x2_node(1))*(x1_node(2)-x1_node(1)) - (x1_node(3)-x1_node(1))*(x2_node(2)-x2_node(1)) );
Na_x1 = ( 0 - (x2_node(3)-x2_node(2)) ) / ( (x2_node(1)-x2_node(2))*(x1_node(3)-x1_node(2)) - (x1_node(1)-x1_node(2))*(x2_node(3)-x2_node(2)) );
Na_x2 = ( (x1_node(3)-x1_node(2)) - 0 ) / ( (x2_node(1)-x2_node(2))*(x1_node(3)-x1_node(2)) - (x1_node(1)-x1_node(2))*(x2_node(3)-x2_node(2)) );
Nb_x1 = ( 0 - (x2_node(1)-x2_node(3)) ) / ( (x2_node(2)-x2_node(3))*(x1_node(1)-x1_node(3)) - (x1_node(2)-x1_node(3))*(x2_node(1)-x2_node(3)) );
Nb_x2 = ( (x1_node(1)-x1_node(3)) - 0 ) / ( (x2_node(2)-x2_node(3))*(x1_node(1)-x1_node(3)) - (x1_node(2)-x1_node(3))*(x2_node(1)-x2_node(3)) );
Nc_x1 = ( 0 - (x2_node(2)-x2_node(1)) ) / ( (x2_node(3)-x2_node(1))*(x1_node(2)-x1_node(1)) - (x1_node(3)-x1_node(1))*(x2_node(2)-x2_node(1)) );
Nc_x2 = ( (x1_node(2)-x1_node(1)) - 0 ) / ( (x2_node(3)-x2_node(1))*(x1_node(2)-x1_node(1)) - (x1_node(3)-x1_node(1))*(x2_node(2)-x2_node(1)) );
u = Na * u_node(1) + Nb * u_node(2) + Nc * u_node(3);
u_x1 = Na_x1 * u_node(1) + Nb_x1 * u_node(2) + Nc_x1 * u_node(3);
u_x2 = Na_x2 * u_node(1) + Nb_x2 * u_node(2) + Nc_x2 * u_node(3);
end
% problem 1.2, L2, L_inf, H1 norm
function [L2,L_inf,H1] = error_norm(coord,connect,connect_number,pressure)
[L2,L_inf,H1] = deal(0,0,0);
for lmn = 1 : connect_number
    % get the position of the three nodes of the triangle element
    x1_node = coord(connect(lmn,1:3),1);
    x2_node = coord(connect(lmn,1:3),2);
    u_node = pressure(connect(lmn,1:3));
    % obtain the edge and area of the elements
    edge_element = [sqrt((x2_node(1)-x2_node(2))^2 + (x1_node(1)-x1_node(2))^2), sqrt((x2_node(1)-x2_node(3))^2 + (x1_node(1)-x1_node(3))^2), sqrt((x2_node(2)-x2_node(3))^2 + (x1_node(2)-x1_node(3))^2)];
    area_element = triangle_area(edge_element);
    % obtain the position of the midpoints
    x1_mid = [(x1_node(1) + x1_node(2))/2 , (x1_node(1) + x1_node(3))/2 , (x1_node(2) + x1_node(3))/2];
    x2_mid = [(x2_node(1) + x2_node(2))/2 , (x2_node(1) + x2_node(3))/2 , (x2_node(2) + x2_node(3))/2];
    for i = 1 : 3 % integration over all integration points
        % the numerical solution of u_mid, u_x1_mid and u_x2_mid obtained
        [u_mid,u_x1_mid,u_x2_mid] = u_predict(x1_mid(i),x2_mid(i),x1_node,x2_node,u_node);
        % the exact solution u_mid_Motz
%        u_mid_Motz = L_shape(x1_mid(i),x2_mid(i));
%        dx = sqrt(2) * 10^(-8); % increment to calculate the exact derivatives
%        u_mid_Motz_x1 = (L_shape(x1_mid(i) - dx, x2_mid(i) + 0 ) - u_mid_Motz)/(-dx);
%        u_mid_Motz_x2 = (L_shape(x1_mid(i) + 0 , x2_mid(i) - dx) - u_mid_Motz)/(-dx);
%
        u_mid_temp = L_exact([x1_mid(i),x2_mid(i)]);
        u_mid_Motz = u_mid_temp(1);
        u_mid_Motz_x1 = u_mid_temp(2);
        u_mid_Motz_x2 = u_mid_temp(3);
%
        u_error = abs(u_mid_Motz - u_mid);
        u_x1_error = abs(u_x1_mid - u_mid_Motz_x1);
        u_x2_error = abs(u_x2_mid - u_mid_Motz_x2);
        % updata L2 norm error
        L2 = L2 + u_error^2/3 * area_element;
        % updata Linf norm error
        if abs(u_error > L_inf)
            L_inf = u_error;
        end
        % updata H1 norm error
        dH1 = (u_error^2 + u_x1_error^2+ u_x2_error^2)/3 * area_element;
        H1 = H1 + (u_error^2 + u_x1_error^2+ u_x2_error^2)/3 * area_element;
    end
end
L2 = sqrt(L2);
H1 = sqrt(H1);
end
% problem 1.1, pointwise error contours
function err_contour(coord,connect,pressure,nnode)
disp("generate error contour...");
pressure_exact = zeros(nnode,1);
pressure_error = zeros(nnode,1);
for i = 1 : nnode
%    pressure_exact(i) = L_shape(coord(i,1),coord(i,2));
%
    temp = L_exact(coord(i,:));
    pressure_exact(i) = temp(1);
%
    pressure_error(i) = abs(pressure_exact(i) - pressure(i));
end
trimesh(connect,coord(:,1),coord(:,2),pressure_error); % hand in this plot
xlabel('x'), ylabel('y'), zlabel('u(x,y)'), az = 70; el = 60; view(az, el);
end
% area of a triangle
function area_element = triangle_area(edge_element)
a = edge_element(1);b = edge_element(2);c = edge_element(3);
area_element = 1/4 * sqrt((a+b+c) * (a+b-c) * (a+c-b) * (b+c-a));
end
%==================================end=====================================