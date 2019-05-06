%% Scientific Computing Project 2D Poisson Eqn.(AP02-2),   Daniel Bosquez
% Gauss Seidel:
clc
close all
fprintf('Running\n') % Message to note code has started/currently running

% CHECKPOINTING is commented, uncomment to enable
% if exist( 'checkpointGS.mat','file' ) % If a checkpoint file exists, load it
%     fprintf('Checkpoint file found - Loading\n');
%     load('checkpointGS.mat')
% 
% else %otherwise, start from the beginning
%     fprintf('No checkpoint file found - starting from beginning\n');
%     mysum=0;
%     countmin=1;
% end

% Define grid

iter = 32000; % Enter # of iterations (32k-34k for best soln.)
N = 600; % Enter # of grid intervals
h = (2*pi)/(N+1); % grid step dx = dy


% Specify Bounds

ax = 0;     % x = 0
ay = ax;    % y = 0
bx = 2*pi;  % x = L
by = bx;    % y = L

st = 1:N+2; % number of x and y steps
len = length(st);
j = st'; % x step vector
k = st; % y step vector
totl = len*len; % Total number of elements
endbc = totl-len+1; 

xval=(h.*j-1); % x values for plotting purposes
yval=(h.*k-1); % y values for plotting purposes


% Initialize and vectorize known conditions

F = sin(pi.*(((h.*(j-1))-ax)./(bx-ax)))*cos((pi/2).*((2.*(((h.*(k-1))-ay)./(by-ay)))+1)); % F matrix of known F's for all x and y's
F = F(:);                                 % Vectorizes F matrix
fa = (h.*(j-1)).*((h.*(j-1))-ax).^2;      % BC equation U(x,y=by)
ga = ((h.*(j-1))-ax).^2.*cos((h.*(j-1))); % BC equation U(x,y=ay)


% Create the U solution vector then populate with known conditions.

U = zeros(len);   % initialize solution array, zero element place holders also act as initial values
U(1:len)= ga;                                                     % U(x,y=ay) Boundary condition
U(endbc:totl)=fa;                                                 % U(x,y=by) Boundary condition
U(len,:) = ga(len)+((((h.*(k-1))-ay)/(bx-ay))*(fa(len)-ga(len))); % U(x=bx,y) Boundary condition
%U2 = U(:); % vectorized array ("turned off")


% Commence Gauss Seidel Matrix Solver

for i= 1:iter %countmin:iter <-- CHECKPOINT STEP % loop solving for every i iteration of GS method until convergence of an adequate set of solutions
   for K = 2:len-1 % Cycling through column entries (Y dimension)
     U(1,K)=(.25*(U(2,K)+U(3,K)+U(2,K-1)+U(2,K+1)))+(.25*h*h*F(1+((K-1)*len))); % "Ghost Node" entries for Neumann condition
        for J = 2:len-1 % Cycling through row entries (X dimension)
        U(J,K)= (.25*(U(J-1,K)+U(J+1,K)+U(J,K-1)+U(J,K+1)))+(.25*h*h*F(J+((K-1)*len))); % Discretized, 4 point Gauss Seidel equation
        
        %U2(J+(K-1)*len) =.25*(U(J-1+(K)*len)+U(J+1+((K)*len))+U(J+((K-1)*len))+U(J+((K)*len))); (vectorized discretization form "turned off", problem with Yk value indexing.)
        end
   end
%  countmin = i+1; %<-- CHECKPOINT STEPS  %If we load this checkpoint, we want to start on the next iteration
%     if mod(count,25)==0 
%         %save checkpoint   
%         fprintf('Saving checkpoint\n');
%         save('checkpointGS.mat');
%     end
end
clc
fprintf('Done') % Signal to user operation is complete


%% Figure Generation for Gauss Seidel Visualizations
% Seperated for speed purposes as it is a non-critical feature
figure;
contour3(xval,yval,U,len,'ShowText','off')
xlabel('0 < X < 2pi')
ylabel('0 < Y < 2pi')
zlabel('U(Xj,Yk)')
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
% zlabel('U3')
% figure;
% plot(1:totl,U(:))
% title('Solution U for every dj,dk')
% xlabel('Step jk')
% ylabel('U(Xj,Yk)')
