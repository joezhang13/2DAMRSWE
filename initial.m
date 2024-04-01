clear;

L = 10;
N = 500;
dx = L/N;
x = ((1:N)-0.5)*dx;
[X, Y] = meshgrid(x);
h0 = 6 + 3*exp(-8*(X-5).^2-5*(Y-5).^2);
figure;
contourf(X,Y,h0,50,'EdgeColor','none');
axis equal;
colorbar;
xlim([0 L]);
ylim([0 L]);
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
set(gca, 'FontSize', 15);

h0x = 6 + 3*exp(-8*(x-5).^2);
h0y = 6 + 3*exp(-5*(x-5).^2);
figure;
plot(x, h0x, 'b-', 'LineWidth', 1.5);
hold on
plot(x, h0y, 'r-', 'LineWidth', 1.5);
xlabel('$x, y$','Interpreter','latex');
ylabel('$h$','Interpreter','latex');
legend('$h(x,y=5)$','$h(x=5,y)$','Interpreter','latex');
set(gca, 'FontSize', 15);