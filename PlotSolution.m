function PlotSolution(X,Y,Z,T)
  figure;
  %pcolor(X,Y,Z); 
  surfl(X,Y,Z); 
  axis([0 1 0 1 0 1]); caxis([0 1]); %shading interp;
  xlabel('x','fontsize',20); ylabel('y','fontsize',20);
  zlabel(strcat('u(x,y,',num2str(T),')'),'fontsize',20);
  title(strcat('t=',num2str(T)),'fontsize',20);
  grid off
return