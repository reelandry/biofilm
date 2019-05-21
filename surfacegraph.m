% creates a titled, single surface plot of each state variable
% then saves the result according to t = T

hf = figure;
set(hf,'position',[10 190 640 500]);
filenameS = strcat('Sr', num2str(paramcond), '_', strrep(num2str(t), '.', '_'), '.eps');
surf(XX,YY,S','FaceColor','interp',...
    'EdgeColor','k',...
    'FaceLighting','phong');
title(strcat('Surface Plot of S at T = ', num2str(t)));%,'Interpreter','LaTex');
xlabel('x');%,'Interpreter','LaTex');
ylabel('y');%,'Interpreter','LaTex');
zlabel(strcat('S(x, y, ', num2str(t), ')'));%,'Interpreter','LaTex');
%daspect([5 5 0.1]);
axis([0, 1, 0, 1, 0, tmaxS])
view(-50,30);
camlight left;
saveas(gcf,filenameS, 'psc2')
closereq

hf = figure;
set(hf,'position',[10 190 640 500]);
filenameX = strcat('Xr', num2str(paramcond), '_', strrep(num2str(t), '.', '_'), '.eps');
surf(XX,YY,X','FaceColor','interp',...
    'EdgeColor','k',...
    'FaceLighting','phong');
title(strcat('Surface Plot of X at T = ', num2str(t)));%,'Interpreter','LaTex');
xlabel('x');%,'Interpreter','LaTex');
ylabel('y');%,'Interpreter','LaTex');
zlabel(strcat('X(x, y, ', num2str(t), ')'));%,'Interpreter','LaTex');
%daspect([5 5 0.1]);
axis([0,1, 0, 1, 0, tmaxX])
view(-50,30);
camlight left;
saveas(gcf,filenameX, 'psc2')
closereq

hf = figure;
set(hf,'position',[10 190 640 500]);
filenameI = strcat('Ir', num2str(paramcond), '_', strrep(num2str(t), '.', '_'), '.eps');
surf(XX,YY,I','FaceColor','interp',...
    'EdgeColor','k',...
    'FaceLighting','phong');
title(strcat('Surface Plot of I at T = ', num2str(t)));%,'Interpreter','LaTex');
xlabel('x');%,'Interpreter','LaTex');
ylabel('y');%,'Interpreter','LaTex');
zlabel(strcat('I(x, y, ', num2str(t), ')'));%,'Interpreter','LaTex');
%daspect([5 5 0.1]);
axis([0,1, 0, 1, 0, tmaxI])
view(-50,30);
camlight left;
saveas(gcf,filenameI, 'psc2')
closereq

hf = figure;
set(hf,'position',[10 190 640 500]);
filenameE = strcat('Er', num2str(paramcond), '_', strrep(num2str(t), '.', '_'), '.eps');
surf(XX,YY,E','FaceColor','interp',...
    'EdgeColor','k',...
    'FaceLighting','phong');
title(strcat('Surface Plot of E at T = ', num2str(t)));%,'Interpreter','LaTex');
xlabel('x');%,'Interpreter','LaTex');
ylabel('y');%,'Interpreter','LaTex');
zlabel(strcat('E(x, y, ', num2str(t), ')'));%,'Interpreter','LaTex');
%daspect([5 5 0.1]);
axis([0,1, 0, 1, 0, tmaxE])
view(-50,30);
camlight left;
saveas(gcf,filenameE, 'psc2')
closereq