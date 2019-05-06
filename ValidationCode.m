% VERIFICATION CODE: Method of Manufactured Solutions
% Daniel Bosquez PSID: 1332758

% **SOR PROCESS STEPS COMMENTED BY DEFAULT**

% if Vx''+Vy''=0 with V(x,y)=sin(x-y) and Vx''+Vy''= -2sin(x-y) = 0 
% Then if Ux''+Uy''= 0 = -2sin(x-y)

% Code will take on the form of SOR method
% To validate GS results, comment SOR terms

clc
clear
close all
fprintf('Running\n') % Message to note code has started/currently running

%addup_checkpoint_rand.m


% Define grid
iter = 32000; % Enter # of iterations for convergence
N = 600; % Enter # of grid intervals (ConvStud: n = 595, w = 1.70524
%w = 1.71; % Enter relaxation coeffient omega (w = 1 is equivelent to Gauss Seidel Solution)
h = (2*pi)/(N+1); % grid step dx = dy

% Specify Bounds
ax = 0;     % x = 0
ay = ax;    % y = 0
bx = 2*pi;  % x = L
by = bx;    % y = L
st = 1:N+2; % number of spatial x and y steps
len = length(st);
j = st'; % x step vector
k = st;  % y step vector
totl = len*len; % Total number of solution entries
endbc = totl-len+1; 

xval=(h.*j-1); % x values for plotting purposes
yval=(h.*k-1); % y values for plotting purposes

%% Manufactured Solution V(x,y)= sin(x-y)
% **SOR Processes commented out by default**
% Initialize and vectorize known conditions
Phi = -2*sin((j-1)-(k-1)); % F matrix of known F's for all x and y's
Phi = Phi(:);                                % Vectorizes F matrix
fa = sin((j-1)-by);     % BC equation U(x,y=by)
ga = sin(j-1);          % BC equation U(x,y=ay)

% Create the U vector & populate with known conditions. Dirichlet data only
Umanu = zeros(len);    % initialize solution array, zero element place holders also act as initial values
Umanu(1:len)= ga;      % U(x,y=ay) Boundary Condition
Umanu(endbc:totl)=fa;  % U(x,y=by) Boundary Condition
Umanu(len,:) = sin(2*pi-(k-1)); % U(bx,y) Boundary Condition
Umanu(1,:) = sin(-(k-1));       % U(ax,y) Boundary Condition
%preU = Umanu;  % initial values for Ujkn-1 (Previous iteration solution)


% Commence SOR Gauss-Seidel Vector solver
for i=1:iter % loop for every i iteration of method until solution convergence
    for K = 2:len-1 % Cycling through column entries (Y dimension)
    % U(1,K)=(.25*(U(2,K)+U(3,K)+U(2,K-1)+U(2,K+1)))+(.25*h*h*F(1+((K-1)*len))); % "Ghost Node" entries for Neumann condition
        for J = 2:len-1 % Cycling through row entries (X dimension)
        Umanu(J,K)= (.25*(Umanu(J-1,K)+Umanu(J+1,K)+Umanu(J,K-1)+Umanu(J,K+1)))+(.25*h*h*Phi(J+((K-1)*len)));     % Explicit Ujkn value for current iterative step n (Gauss Seidel soln)
%         U(J,K)=w*U(J,K)+(1-w)*preU(J,K);    % SOR expression: Implicit Ujkn+1 = w*(Explicit Ujkn)+(1-w)*(Previous Ujkn-1 from last iteration)
%         preU(J,K) = U(J,K);                 % Ujkn-1 term for next n iteration
         
        %U2(J+(K-1)*len) = .25*(U(J-1+(K)*len)+U(J+1+((K)*len))+U(J+((K-1)*len))+U(J+((K)*len))); %(Vectorized discretization form "turned off", problem with Yk value indexing.)
        end
    end   
end
clc
fprintf('Done\n') % Signal to user operation is complete
% Error Check for Final Validaion
Err1 = sum(abs(Umanu-(sin(xval'-yval))))*(1/(J*K));
Errinf = max(abs(Umanu-(sin(xval'-yval))));
AverageAbsErr = mean(Err1)
AverageInfErr = mean(Errinf)

