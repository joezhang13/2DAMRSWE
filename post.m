clear;
output = load('output2.dat');

L = 10;
N = 100;
dxc = L/N;
dxf = L/(2*N);
xc = ((1:N)-0.5)*dxc;
xf = ((1:2*N)-0.5)*dxf;
I = output(:, 1);
idx = output(:, 2);
h = output(:, 3);
id = output(:, 6);

nn = length(I);
xx = zeros(4, nn);
yy = zeros(4, nn);
for i = 1 : nn
    j = idx(i);
    if I(i)
        jc = floor(j/4);
        jf = mod(j, 4);
        xci = xc(mod(jc,N)+1);
        yci = xc(floor(jc/N)+1);
        if mod(jf, 2)
            xi = xci + 0.5*dxf;
        else
            xi = xci - 0.5*dxf;
        end
        if jf > 1
            yi = yci + 0.5*dxf;
        else
            yi = yci - 0.5*dxf;
        end
        xx(:, i) = [-0.5;0.5;0.5;-0.5]*dxf + xi;
        yy(:, i) = [-0.5;-0.5;0.5;0.5]*dxf + yi;
    else
        xi = xc(mod(j,N)+1);
        yi = xc(floor(j/N)+1);
        xx(:, i) = [-0.5;0.5;0.5;-0.5]*dxc + xi;
        yy(:, i) = [-0.5;-0.5;0.5;0.5]*dxc + yi;
    end
end

figure;
% patch(xx, yy, h, 'EdgeColor', 'none');
patch(xx, yy, h);
axis equal;
colorbar;
caxis([5.85 6.25]);
xlim([0 L]);
ylim([0 L]);
xlabel('x');
ylabel('y');
set(gca, 'FontSize', 15);

% figure;
% patch(xx, yy, id);
% axis equal;
% xlim([0 L]);
% ylim([0 L]);
% xlabel('x');
% ylabel('y');
% set(gca, 'FontSize', 15);