c = [1 -2 1 -1];  
x = linspace(-2,4);  
y = c(1)*x.^3+c(2)*x.^2+c(3)*x+c(4) + randn(1,100);
plot(x,y,'.b-')

hold on
c = polyfit(x,y,3);  
yhat = c(1)*x.^3+c(2)*x.^2+c(3)*x+c(4);
plot(x,yhat,'r','linewidth',2)

x0 = 1;
y0 = 10;

x = x(:); %reshape the data into a column vector
y = y(:);
% 'C' is the Vandermonde matrix for 'x'
n = 3; % Degree of polynomial to fit
V(:,n+1) = ones(length(x),1,class(x));
for j = n:-1:1
     V(:,j) = x.*V(:,j+1);
end
C = V;
% 'd' is the vector of target values, 'y'.
d = y;
% There are no inequality constraints in this case, i.e., 
A = [];
b = [];
% We use linear equality constraints to force the curve to hit the required point. In
% this case, 'Aeq' is the Vandermoonde matrix for 'x0'
Aeq = x0.^(n:-1:0);
% and 'beq' is the value the curve should take at that point
beq = y0;
p = lsqlin( C, d, A, b, Aeq, beq )
% We can then use POLYVAL to evaluate the fitted curve
yhat = polyval( p, x );
% Plot original data
plot(x,y,'.b-') 
hold on
% Plot point to go through
plot(x0,y0,'gx','linewidth',4) 
% Plot fitted data
plot(x,yhat,'r','linewidth',2) 
hold off