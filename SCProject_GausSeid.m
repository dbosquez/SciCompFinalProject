%% Scientific Computing Project 2D Poisson Eqn.(AP02-2),   Daniel Bosquez
% Gauss Seidel:
clc
close all



% Define grid

iter = 500; %iterations for convergence
N = 1000; %grid intervals
h = (2*pi)/(N+1); %grid step dx = dy
ax = 0;
ay = ax;
bx = 2*pi;
by = bx;
st = [1:N+2]; % number of x and y steps
len = length(st);
totl = len*len;
endbc = totl-len+1;
j = st';
k = st;
%j=2;
%k=2;
xval=(h.*j-1);
yval=(h.*k-1);
% Initialize and vectorize knowns

F = sin(pi.*(((h.*(j-1))-ax)./(bx-ax)))*cos((pi/2).*((2.*(((h.*(k-1))-ay)./(by-ay)))+1)); % F matrix of known F's for all x and y's
F = F(:); % Vectorizes F matrix
fa = (h.*(j-1)).*((h.*(j-1))-ax).^2;
ga = ((h.*(j-1))-ax).^2.*cos((h.*(j-1)));

%create the U vector then populate with known conditions.

U = zeros(len); % initialize solution array
U(1:len)= ga; % U(x,y=ay) BC
U(endbc:totl)=fa; % U(x,y=by) BC
U(len,:) = ga(len)+((((h.*(k-1))-ay)/(bx-ay))*(fa(len)-ga(len))); % U(bx,y) BC
U2 = U(:); % vectorized array
%U=U(:);
% Commence Gauss Seidel Vector solver

for i=1:iter
 for K = 2:len-1
     for J = 2:len-1
         U(J,K)= (.25*(U(J-1,K)+U(J+1,K)+U(J,K-1)+U(J,K+1)))+(.25*h*h*F(J+((K-1)*len)));
     %U2(J+(K-1)*len) = .25*(U(J-1+(K)*len)+U(J+1+((K)*len))+U(J+((K-1)*len))+U(J+((K)*len)));
     end
 end

end

figure;
contour3(xval,yval,U,len,'ShowText','off')
xlabel('0 < X < 2pi')
ylabel('0 < Y < 2pi')
zlabel('U')


figure;
contour3(144,144,U144,len,'ShowText','off')
xlabel('0 < X < 2pi')
ylabel('0 < Y < 2pi')
zlabel('U')


