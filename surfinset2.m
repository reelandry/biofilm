hf = figure;
set(hf, 'Visible','off');
set(hf, 'position',[0 0 1300 1300]);
set(hf, 'MenuBar', 'none');
set(hf, 'ToolBar', 'none');
g = area(x, St, 'FaceColor', [0 1 0] );
axis([0, 1, 0, 1]);
% Make area transparent
drawnow; pause(0.05);  % for transparency to work
g.Face.ColorType = 'truecoloralpha';
g.Face.ColorData(4) = 255 * 0.4;
axis([0, 1, 0, 1]); hold on;
g = area(x, Xt, 'FaceColor', [0 0 1] );

% Make area transparent
drawnow; pause(0.05);  % for transparency to work
g.Face.ColorType = 'truecoloralpha';
g.Face.ColorData(4) = 255 * 0.3;
axis([0, 1, 0, 1]);
g = area(x, It, 'FaceColor', [1 0 0] );

% Make area transparent
drawnow; pause(0.05);  % for transparency to work
g.Face.ColorType = 'truecoloralpha';
g.Face.ColorData(4) = 255 * 0.2;
axis([0, 1, 0, 1]);
g = area(x, Et, 'FaceColor', [1 0 1] );

% Make area transparent
drawnow; pause(0.05);  % for transparency to work
g.Face.ColorType = 'truecoloralpha';
g.Face.ColorData(4) = 255 * 0.1;
axis([0, 1, 0, 1]);legend('S', 'X', 'I', 'E');
title(['T = ', num2str(t), ' ', 'Neu = ', num2str(Neu)]);
axis([0, 1, 0, 1]);
%hold off;
