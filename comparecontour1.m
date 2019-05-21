function E = comparecontour1(XX, YY, UX, VX, US, VS, M, N)
	
	figure;
	subplot(2, 1, 1);
	contourf(XX, YY, log(abs((VX'-UX')./VX'/eps)));colorbar;
	
	subplot(2, 1, 2);
	contourf(XX, YY, log(abs((VS'-US')./VS'/eps)));colorbar;
	
end