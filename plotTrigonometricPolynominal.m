function plotTrigonometricPolynominal(y)
if ~isvector(y)
    error('>>>Error: input must be a vector.');
end
if iscolumn(y)
    y = y';
end
L = 1;
m = length(y);
deltaL = L / m;
isOddGrid = mod(m,2);
x = (deltaL/2):deltaL:(L-deltaL/2);
n = floor((m+1)/2);
z = fft(y)/m;

a0 = z(1); 
an = 2*real(z(2:n));
bn = -2*imag(z(2:n));
if ~isOddGrid
    anp1 = z(n+1);
end

px = 0:(L/1000):L;
px_biase = px - deltaL/2;
k = 1:length(an);
py = a0 + an*cos(2*pi*k'*px_biase/L) ...
        + bn*sin(2*pi*k'*px_biase/L);
if ~isOddGrid
    py = py + anp1 * cos(2*pi*n*px_biase/L);
end

figure
hold on
xlim([0 L]);
plot(x,y,'ko','LineWidth',1.5);
plot(px,py,'k-','LineWidth',1.5);
box on
set(gca,'FontName','Times New Roman','FontSize',16,'LineWidth',1.5);
hold off

end