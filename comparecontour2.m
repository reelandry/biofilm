function E = comparecontour2(x, UX, VX, US, VS, UI, VI, UE, VE)


%subplot(2, 2, 1);
%plot(x, log(abs((VS'-US')./VS')));
%subplot(2, 2, 2);
semilogx(x, log(abs((VX'-UX')./VX')));
%subplot(2, 2, 3);
%plot(x, log(abs((VI'-UI')./VI')));
%subplot(2, 2, 4);
%plot(x, log(abs((VE'-UE')./VE')));
    hold on;
end