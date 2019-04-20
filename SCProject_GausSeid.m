%% Scientific Computing Project 2D Poisson Eqn.(AP02-2),   Daniel Bosquez
% Gauss Seidel:
clc
close
clear

% Define grid

iter = 5; %iterations for convergence
N = 10; %grid intervals
h = (2*pi)/(N+1); %grid step dx = dy

x = [1:N+1] % number of x steps
y = x;      % number of y steps







