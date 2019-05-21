aviObject = VideoWriter('myVideo.avi');
aviObject.FrameRate = 20; % Optional
open(aviObject);
for i=0:100
    x=cos(0:0.02:10*pi)+(0:0.02:10*pi)/pi*(1+sin(i/100*pi*4))*2;
    y=sin((0:0.02:10*pi));
    plot(x,y,'linewidth',2);
    axis([-3,50,-1.2,1.2]);
    drawnow;
    set(gcf,'Visible', 'off', 'color', [1,1,1], 'Position', [1, 1, 1024,768]); 
    I=export_fig('-nocrop');
    F = im2frame(I);
	            frame = getframe(gcf);
            writeVideo(aviObject, frame);

 %   aviObject = addframe(aviObject,F);
    i
end
close;
aviObject = close(aviObject);