%% Scientific Computing Project 2D Poisson Eqn.(AP02-2),   Daniel Bosquez
% Gauss Seidel:
clc
close all

%addup_checkpoint_rand.m


%Checkpoint every 25 iterations


% Define grid

iter = 6000; %iterations for convergence
N = 600; %grid intervals
h = (2*pi)/(N+1); %grid step dx = dy
w = 1; % relaxation coeffient omega (w = 1 is equivelent to Gauss Seidel Solution)

% Specify Bounds
ax = 0;
ay = ax;
bx = 2*pi;
by = bx;
st = 1:N+2; % number of x and y steps
len = length(st);
j = st'; % x step vector
k = st;  % y step vector
totl = len*len;
endbc = totl-len+1;
xval=(h.*j-1); % x values for plotting purposes
yval=(h.*k-1); % y values for plotting purposes

% Initialize and vectorize knowns

fprintf('Running\n') % Message to note code has started/currently running
F = sin(pi.*(((h.*(j-1))-ax)./(bx-ax)))*cos((pi/2).*((2.*(((h.*(k-1))-ay)./(by-ay)))+1)); % F matrix of known F's for all x and y's
F = F(:);                               % Vectorizes F matrix
fa = (h.*(j-1)).*((h.*(j-1))-ax).^2;    % BC equation U(x,y=by)
ga = ((h.*(j-1))-ax).^2.*cos((h.*(j-1)));%BC equation U(x,y=ay)

% Create the U vector then populate with known conditions.

U = zeros(len);    % initialize solution array, zero element place holders also act as initial values
U(1:len)= ga;      % U(x,y=ay) BC
U(endbc:totl)=fa;  % U(x,y=by) BC
U(len,:) = ga(len)+((((h.*(k-1))-ay)/(bx-ay))*(fa(len)-ga(len))); % U(bx,y) BC
preU = U;  % initial value for Ujkn-1 (Previous iteration solution)
%U2 = U(:); % vectorized array ("turned off")


% Commence SOR Gauss-Seidel Vector solver

for i=1:iter % loop for every i iteration of method until solution convergence
    for K = 2:len-1 % Cycling through column entries (Y dimension)
    U(1,K)=(.25*(U(2,K)+U(3,K)+U(2,K-1)+U(2,K+1)))+(.25*h*h*F(1+((K-1)*len))); % "Ghost Node" entries for Neumann condition
        for J = 2:len-1 % Cycling through row entries (X dimension)
        U(J,K)= (.25*(U(J-1,K)+U(J+1,K)+U(J,K-1)+U(J,K+1)))+(.25*h*h*F(J+((K-1)*len))); % Explicit Ujkn value for current iterative step n (Gauss Seidel soln)
        U(J,K)=w*U(J,K)+(1-w)*preU(J,K); % SOR expression: Implicit Ujkn+1 = w*(Explicit Ujkn)+(1-w)*(Previous Ujkn-1 from last iteration)
        preU(J,K) = U(J,K); % Ujkn-1 term for next n iteration
         
        %U2(J+(K-1)*len) = .25*(U(J-1+(K)*len)+U(J+1+((K)*len))+U(J+((K-1)*len))+U(J+((K)*len))); %(Vectorized discretization form "turned off", problem with Yk value indexing.)
        end
    end   
end
clc
fprintf('Done\n') % Signal to user operation is complete


%% Figure Generation for SOR Visualizations
% Seperated for speed purposes as it is a non-critical feature

% Contour plot to visualize where solutions for U(x,y) are defined. Plot
% aids in determining the level of resolution for a given set of iterative
% parameters when comparing different sets.

 figure;
 contour3(xval,yval,U,len,'ShowText','off')
 xlabel('0 < X < 2pi')
 ylabel('0 < Y < 2pi')
 zlabel('U(Xj,Yk)')

 % Comparision Plots to compare against other iterations and steps (Used to
 % determine grid convergence)
 
%figure;
% contour3(x500,y500,U500,len,'ShowText','off')
% xlabel('0 < X < 2pi')
% ylabel('0 < Y < 2pi')
% zlabel('U1')
%  figure;
%  contour3(x200,y200,U200,len,'ShowText','off')
%  xlabel('0 < X < 2pi')
%  ylabel('0 < Y < 2pi')
%  zlabel('U2')
%figure;
% contour3(x62,y62,UN62,len,'ShowText','off')
% xlabel('0 < X < 2pi')
% ylabel('0 < Y < 2pi')
% % zlabel('U3')

%2D plot (Unused)

%  figure;
%  plot(1:totl,U(:))
%  title('Solution U for every dj,dk')
%  xlabel('Step jk')
%  ylabel('U(Xj,Yk)')
