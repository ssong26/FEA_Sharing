% =====================29-Apr-2018 Final Project============================
% Siyuan_Song Final Project for APMA 2560
% 29-Apr-Problem 1
%
% Program Description------------------------------------------------------
%
% 1. This program is to solve the Motz problem for a cube.
%
% 2. The following functions are included in the program
% 2.1 Pde_Final_1()-----------------
%     The main function
% 2.2 grid_generation()-------------
%     The uniform grid is generated, with the (x_number,y_number) as input
%     [coord,connect,connect_number,nnode] as the output
% 2.3 
%
% 3. Parameter introduction
% coeff------------------------------wave number of Helmoholz equation
% x_min------------------------------minimum value of x of the structure
% y_min------------------------------minimum value of y of the structure
% x_max------------------------------maximum value of x of the structure
% y_max------------------------------maximum value of y of the structure
% u_1  ------------------------------value of u at the direchilit-boundary 1
% u_3  ------------------------------value of u at the direchilit-boundary 1
% x_number --------------------------number of node per unit length in x_direction
% y_number --------------------------number of node per unit length in y_direction
% coord------------------------------(:,2)matrix, coord(:,1) record the x
% position, coord(:,2) record the y position of the nodes of the elements
% connect----------------------------(:,3)matrix, record the series number
% of the three nodes of each trianglular elements.
% connect_number---------------------number of elements in the system
% nnode------------------------------number of nodes in the system
% pressure---------------------------(:,1)vector, fem solutions of the aim function u
% 
%
% ==========================================================================
% =========================The Main Function================================
%
function PDE_Final_1()
clear;clc;close all
[xy_min,xy_max,xy_number] = deal(10,100,10);
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
    [coeff,x_min,x_max,y_min,y_max,u_1,u_3] = deal(0,-1,1,0,1,0,500);
    %
    %----------------------------Part 2 Generating Grid------------------------
    %
    % ---Grid parameters
    [x_number,y_number,alpha,layer] = deal(xy_vector(i),xy_vector(i),1,floor(xy_vector(i)*4/5));
    % ---Obtain the grid
    % [coord,connect,connect_number,nnode,rh] = grid_iteractive(x_min,x_max,y_min,y_max,xy_vector(i),alpha,layer);
    [coord,connect,connect_number,nnode,rh] = grid_generation(x_min,x_max,y_min,y_max,x_number,y_number,alpha,layer);
    %
    % ---------------------Part 3 Finite element analysis-----------------------
    %
    pressure = fem_Motz(coord,connect,connect_number,nnode,u_1,u_3,coeff);
    %
    % ========================Compare the result================================
    %
    % ================================
    % the contour of FEM
    % fem_contour(x_min,x_max,y_min,y_max,coord,pressure);
    % ================================
    % the contour of exact
    % exact_contour(x_min,x_max,y_min,y_max,coord,nnode);
    % ===============================
    % comparison figure of FEM and theory at the line y = y_line
    % y_yline = 0.5;
    % compare_yline(x_min,x_max,y_min,y_max,x_number,y_number,coord,pressure,y_yline);
    %
    % program 1, the pointwise error contour
    % err_contour(x_min,x_max,y_min,y_max,coord,pressure,nnode);
    % derr_contour(x_min,x_max,y_min,y_max,coord,connect,connect_number,pressure);
    %
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
% figure;
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
%
% grid generation function. (wrong)
function [coord,connect,connect_number,nnode,rh] = grid_generation(x_min,x_max,y_min,y_max,x_number,y_number,alpha,layer)
nnode = (2 * x_number - 1) * y_number;    % number of node in the structure
coord = zeros(nnode,2);                   % coordinates of the node

if alpha == 1
    x_temp = linspace(x_min,x_max,(2 * x_number - 1));     % record the x position for nodes
    y_temp = linspace(y_min,y_max,y_number);     % record the y position for nodes
    % ---Generate the node and the basis connect
    rh = (y_max-y_min)/(y_number-1);
else
    if layer == 0  %  
        hx = (alpha - 1)/(alpha^(x_number-1)-1);
        hy = (alpha - 1)/(alpha^(y_number-1)-1);
        x_temp = zeros(1,(2*x_number-1));
        y_temp = zeros(1,y_number);
        for i = 2 : x_number
            x_temp(x_number - i + 1) = x_temp(x_number - i + 2) - alpha^(i-2)*hx;
            x_temp(x_number + i - 1) = x_temp(x_number + i - 2) + alpha^(i-2)*hx;
        end
        for i = 2 : y_number
            y_temp(i) = y_temp(i-1) + alpha^(i-2)*hy;
        end
    else
        hx = 1/((alpha^layer-1)/(alpha-1)+(x_number-1-layer)*alpha^(layer-1));
        hy = 1/((alpha^layer-1)/(alpha-1)+(y_number-1-layer)*alpha^(layer-1));
        x_temp = zeros(1,(2*x_number-1));
        y_temp = zeros(1,y_number);
        for i = 2 : (layer+1)
            x_temp(x_number - i + 1) = x_temp(x_number - i + 2) - alpha^(i-2)*hx;
            x_temp(x_number + i - 1) = x_temp(x_number + i - 2) + alpha^(i-2)*hx;
        end 
        for i = (layer+2) : x_number
            x_temp(x_number - i + 1) = x_temp(x_number - i + 2) - alpha^(layer-1)*hx;
            x_temp(x_number + i - 1) = x_temp(x_number + i - 2) + alpha^(layer-1)*hx;
        end
        for i = 2 : (layer+1)
            y_temp(i) = y_temp(i - 1) + alpha^(i-2)*hy;
        end 
        for i = (layer+2) : y_number
            y_temp(i) = y_temp(i - 1) + alpha^(layer-1)*hy;
        end
        rh = alpha^(layer-1)*hx;
    end
end
% generate the grid
for i = 1 : (2 * x_number - 1)
    for j = 1 : y_number
        node_count = (i - 1) * y_number + j;
        coord(node_count,1)= x_temp(i);
        coord(node_count,2)= y_temp(j);
    end
end
connect=delaunay(coord(:,1),coord(:,2));%generate the connectivity array
%
connect_number=size(connect,1);
%
% figure;
% triplot(connect,coord(:,1),coord(:,2),'r');
end
% grid generatoin iteractive function. (wrong)
function [coord,connect,connect_number,nnode,rh] = grid_iteractive(x_min,x_max,y_min,y_max,r_number,alpha,layer)
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
    coord(number,1) = - r_vector(i+1);
    coord(number,2) = 0;
    number = number + 1;
    r_temp = linspace(increment,r_vector(i+1),temp_number);
    for j = 1 : temp_number
        coord(number,1) = - r_vector(i+1);
        coord(number,2) = r_temp(j);
        number = number + 1;
    end
    for j = 1 : temp_number
        coord(number,1) = - r_vector(i+1) + r_temp(j);
        coord(number,2) = r_vector(i+1);
        number = number + 1;
    end
    for j = 1 : temp_number
        coord(number,1) = 0 + r_temp(j);
        coord(number,2) = r_vector(i+1);
        number = number + 1;
    end
    for j = 1 : temp_number
        coord(number,1) = r_vector(i+1);
        coord(number,2) = r_vector(i+1) - r_temp(j);
        number = number + 1;
    end
end
%
nnode = length(coord(:,1));
connect=delaunayTriangulation(coord(:,1),coord(:,2));%generate the connectivity array
%
connect_number=size(connect,1);
%
% figure;
% triplot(connect,coord(:,1),coord(:,2),'r');
rh = rh * alpha^(layer-1);
end
% finite element program
function pressure = fem_Motz(coord,connect,connect_number,nnode,u_1,u_3,coeff)
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
%{
figure
spy(Stif)
p = symrcm(Stif); 
figure
spy(Stif(p,p))
%}
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
fixnodes = [];
temp_number = 1;
for j = 1 : nnode
    % Curve_1 u = u_1
    if coord(j,1) <= 0 && abs(coord(j,2))<10^(-10)
        fixnodes(temp_number,1) = j;
        fixnodes(temp_number,2) = u_1;
        temp_number = temp_number + 1;
    end
    if abs((coord(j,1) - 1))<10^(-10)
        fixnodes(temp_number,1) = j;
        fixnodes(temp_number,2) = u_3;
        temp_number = temp_number + 1;
    end
end
fixnodes_number = length(fixnodes(:,1));
% figure
% scatter(coord(fixnodes(:,1),1),coord(fixnodes(:,1),2));
%
% ==================== Assemble Boundary condition ========================
%
% Modify the global stiffness and residual to include constraints
% The symmetry must be satisfied
%
resid = zeros(nnode,1);
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
Stif_sparse = sparse(Stif);
pressure=Stif_sparse\resid;
end
% contour of the fem solutions
function fem_contour(x_min,x_max,y_min,y_max,coord,pressure)
x_mesh_number = 100;
y_mesh_number = 100;
xq = linspace(x_min,x_max,x_mesh_number);yq = linspace(y_min,y_max,y_mesh_number);
[Xq,Yq] = meshgrid(xq,yq);
Zq=griddata(coord(:,1),coord(:,2),pressure,Xq,Yq,'v4');
figure;
gca = pcolor(Xq,Yq,Zq);
set(gca,'LineStyle','none')
colorbar;
end
% contour of the exact solutions
function exact_contour(x_min,x_max,y_min,y_max,coord,nnode)
x_mesh_number = 100;
y_mesh_number = 100;
xq = linspace(x_min,x_max,x_mesh_number);yq = linspace(y_min,y_max,y_mesh_number);
[Xq,Yq] = meshgrid(xq,yq);
pressure_exact = zeros(nnode,1);
for i = 1 : nnode
    pressure_exact(i) = Motz(coord(i,1),coord(i,2));
end
Zq=griddata(coord(:,1),coord(:,2),pressure_exact,Xq,Yq,'v4');
figure;
gca = pcolor(Xq,Yq,Zq);
colorbar;
set(gca,'LineStyle','none');
end
% compare the value of u on the position y = y_line
function compare_yline(x_min,x_max,y_min,y_max,x_number,y_number,coord,pressure,y_yline)
figure;
% theory
x_yline_number = 200;
x_yline = linspace(x_min,x_max,x_yline_number);
u_yline_exact = zeros(1,x_yline_number);
for i = 1 : x_yline_number
    u_yline_exact(i) = Motz(x_yline(i),y_yline);
end
plot(x_yline,u_yline_exact,'LineWidth',4);hold on
% fem
y_position = floor(abs((y_yline-0.000000001))/(y_max-y_min)*y_number) + 1; 
scatter(coord(y_position:y_number:(2*x_number-1)*y_number,1),pressure(y_position:y_number:(2*x_number-1)*y_number),150)
legend('theoretical prediction','numerical simulation')
xlabel('X-direction');ylabel('u(x,y)')
end
% problem 1.1, pointwise error contours
function err_contour(x_min,x_max,y_min,y_max,coord,pressure,nnode)
pressure_exact = zeros(nnode,1);
pressure_error = zeros(nnode,1);
for i = 1 : nnode
    pressure_exact(i) = Motz(coord(i,1),coord(i,2));
    pressure_error(i) = abs(pressure_exact(i) - pressure(i));
end
x_mesh_number = 100;
y_mesh_number = 100;
xq = linspace(x_min,x_max,x_mesh_number);yq = linspace(y_min,y_max,y_mesh_number);
[Xq,Yq] = meshgrid(xq,yq);
Zq=griddata(coord(:,1),coord(:,2),pressure_error,Xq,Yq,'v4');
figure
gca=pcolor(Xq,Yq,Zq);
xlabel('x');ylabel('y');
h=colorbar;caxis([1 5]);
h.Ticks=1:5;
h.TickLabels = [10^(-3),10^(-2),10^(-1),10^(0),10^(1)];
set(gca,'LineStyle','none')
end
function derr_contour(x_min,x_max,y_min,y_max,coord,connect,connect_number,pressure)
ux_error = zeros(connect_number,1);
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
    dx = 10^(-6);
    u_center_Motz = Motz(x1_center,x2_center);
    u_center_Motz_x1 = (Motz(x1_center + dx, x2_center + 0 ) - u_center_Motz)/dx;
    u_center_Motz_x2 = (Motz(x1_center + 0 , x2_center + dx) - u_center_Motz)/dx;
    % 
    ux_error(lmn) = sqrt( (u_x1- u_center_Motz_x1)^2 + (u_x2-u_center_Motz_x2)^2 );
    % ux_error(lmn) = sqrt( u_center_Motz_x1^2 + u_center_Motz_x2^2 );
    ux_coord(lmn,1) = x1_center;
    ux_coord(lmn,2) = x2_center;
end
x_mesh_number = 100;
y_mesh_number = 100;
xq = linspace(x_min,x_max,x_mesh_number);yq = linspace(y_min,y_max,y_mesh_number);
[Xq,Yq] = meshgrid(xq,yq);
Zq=griddata(ux_coord(:,1),ux_coord(:,2),ux_error,Xq,Yq,'v4');
figure
gca=pcolor(Xq,Yq,Zq);
xlabel('x');ylabel('y');
h=colorbar;caxis([1 5]);
h.Ticks=1:5;
h.TickLabels = [10^(-3),10^(-2),10^(-1),10^(0),10^(1)];
set(gca,'LineStyle','none')
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
        u_mid_Motz = Motz(x1_mid(i),x2_mid(i));
        dx = 10^(-10); % increment to calculate the exact derivatives
        u_mid_Motz_x1 = (Motz(x1_mid(i) + dx, x2_mid(i) + 0 ) - u_mid_Motz)/dx;
        u_mid_Motz_x2 = (Motz(x1_mid(i) + 0 , x2_mid(i) + dx) - u_mid_Motz)/dx;
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
        H1 = H1 + (u_error^2 + u_x1_error^2+ u_x2_error^2)/3 * area_element;
    end
end
L2 = sqrt(L2);
H1 = sqrt(H1);
end
% area of a triangle
function area_element = triangle_area(edge_element)
a = edge_element(1);b = edge_element(2);c = edge_element(3);
area_element = 1/4 * sqrt((a+b+c) * (a+b-c) * (a+c-b) * (b+c-a));
end
% function determine the value of u at the point (x1,x2)
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
    kel = kel - 1*(C-coeff^2*T);
end
end
% Exact solution function, we can obtain the exact solution at the point (x,y)
function s = Motz(x,y)
% initially, determine the coefficient
coeff = [401.1624537452,87.6559201951,17.2379150794,-8.0712152597,1.4402727170,0.3310548859,0.2754373445,-0.0869329945,0.0336048784,0.0153843745,0.0073023017,-0.0031841139,0.0012206461,0.0005309655,0.0002715122,-0.0001200463,0.0000505400,0.000023167,0.000011535,-0.000005295];
%first, obtain the radial coordinates
r = sqrt(x^2+y^2);
if x>0
    theta = atan(y/x);
end
if x<0
    theta = pi - atan(-y/x);
end
if x==0
    theta = pi/2;
end
% second, obtain the real solution
s = 0;
for k = 1 : 20
    s = s + coeff(k) * (r^((2*k-1)/2)) * cos((2*k-1)*theta/2);
end
end
%==================================end=====================================