%% Scientific Computing Project 2D Poisson Eqn.(AP02-2),   Daniel Bosquez
% Gauss Seidel:
clc
close
clear

% Define grid

iter = 5; %iterations for convergence
N = 10; %grid intervals
h = (2*pi)/(N+1); %grid step dx = dy
ax = 0;
ay = ax;
bx = 2*pi;
by = bx;
st = [1:N+2] % number of x and y steps
j = st';
k = st;
%j=2;
%k=2;

% Initialize and vectorize knowns

F = sin(pi.*(((h.*(j-1))-ax)./(bx-ax)))*cos((pi/2).*((2.*(((h.*(k-1))-ay)./(by-ay)))+1)); % F matrix of known F's for all x and y's
F = F(:); % Vectorizes F matrix
fa = (h.*(j-1)).*((h.*(j-1))-ax).^2;
ga = ((h.*(j-1))-ax).^2.*cos((h.*(j-1)));

% Next up I need to create the U vector then populate with known
% conditions.




