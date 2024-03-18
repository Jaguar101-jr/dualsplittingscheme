
%% Codes by Hussein A. H. Muhammed @School of Stat. & Math. NWPU, Xi'an, China March 2024©.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Poisson Equ. solution A*p=b in dual-splitting scheme

clear;
clc;

Nx = 100; % number of points in the x-direction
Ny = 100; % number of points in the y-direction
Lx = 1.0; % domain length in x
Ly = 1.0; % domain length in y

dx = Lx / (Nx - 1); % grid spacing in x
dy = Ly / (Ny - 1); % grid spacing in y

assert(dx == dy, 'Grid spacings dx and dy should be equal for this example.');
S = zeros(Nx, Ny);


for i = 1:Nx
    for j = 1:Ny
        S(i,j) = sin(pi * (i-1) * dx) * cos(pi * (j-1) * dy);
    end
end
% Number of unknowns
N = Nx * Ny;

% Initialize the matrix A and vector b
A = sparse(N, N); % Sparse matrix for efficiency
b = zeros(N, 1);

% Fill A and b
for i = 1:Nx
    for j = 1:Ny
        
        index = (j-1) * Nx + i;
        
      
        b(index) = -dx^2 * S(i,j);
        
        % Diagonal element
        A(index, index) = -4;
        
        % Off-diagonal elements
        if i > 1
            A(index, index - 1) = 1; 
        end
        if i < Nx
            A(index, index + 1) = 1; 
        end
        if j > 1
            A(index, index - Nx) = 1; 
        end
        if j < Ny
            A(index, index + Nx) = 1; 
        end
    end
end
% Solve the linear system
p = A\b;

% Reshape p back into a 2D grid
P = reshape(p, [Nx, Ny]);
% Create a meshgrid for plotting
[x, y] = meshgrid(0:dx:Lx, 0:dy:Ly);

% Plot the pressure field
surf(x, y, P');
xlabel('x');
ylabel('y');
zlabel('Pressure');
title('Pressure Field');

