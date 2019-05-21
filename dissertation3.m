%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Biofilm Experiment 3                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning off; %#ok<*WNOFF>
clc;
clear all;

% begin recording time of total program
tic;

% end time
T = 1;

% evaluate using predefined paramaters and initial conditions
p = 3;

% display total amount of memory used
checkmemory = 0;

% flag to create and save animation and not save images
animate = 1;

% flag to create a window showing estimated time left
waitbarflag = 0;

% loop through two conditions only changing one setting
for loop = 2:2
	
	% establish intial spatial and temporal resolution
	dx = 0.05;
	dy = 0.05;
	dt = 0.01;
	
	% this performs loop over Neumann and Dirichlet conditions
	if(loop == 1)
		
		Neu = 0;
		
	else
		
		Neu = 1;
		
	end
	
	% mesh grid settings, number of node points
	Rx = dt/dx/dx;
	Ry = dt/dy/dy;
	x = 0:dx:1;
	y = 0:dy:1;
	t = 0:dt:T;
	M = length(x)-1;
	N = length(y)-1;
	[XX, YY] = meshgrid(x, y);
	
	% predefined settings for parameters
	switch p
		
		case 1
			
			% set 1:
			mu = 1;
			yh = 0.35;
			kl = 0.01;
			ks = 0.8;
			ke = 0.03;
			ki = 0.4;
			ye = 0.03;
			alpha = 1;
			beta = 6;
			d1 = 0.0002;
			d2 = 0.0001;
			timevec = [0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.0, 1.5, 2, ...
				2.5, 3, 4, 5, 6, 7, 8, 9, 10];
			
		case 2
			
			% set 2:
			mu = 0.02;
			yh = 0.07;
			kl = 0.1;
			ks = 0.2;
			ke = 0.01;
			ki = 2;
			ye = 0.8;
			alpha = 2;
			beta = 2;
			d1 = 0.0002;
			d2 = 0.0001;
			timevec = [0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.0];
			
		case 3
			
			% set 3:
			mu = 10;
			yh = 0.5;%0.05
			ye = 0;
			ks = 0.2;
			ki = 2;
			ke = 0;
			kl = 0.1;
			alpha = 2;
			beta = 2;
			d1 = 0.002;
			d2 = 0.001;
			timevec = [0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.0, 1.5, 2, ...
				2.5, 3, 4, 5, 6, 7, 8, 9, 10];
		
		case 4
			
			% set 4:
			mu = 1;
			yh = 0.35;
			kl = 0.01;
			ks = 0.8;
			ke = 0.03;
			ki = 0.4;
			ye = 0.03;
			alpha = 2;
			beta = 2;
			d1 = 1;
			d2 = 1;
			timevec = [0, 0.5, 1, 1.5, 2, 2.5];
			
		case 5
			
			% set 5:
			mu = 0.02;
			yh = 0.07;
			kl = 0.1;
			ks = 0.2;
			ke = 0.01;
			ki = 2;
			ye = 0.8;
			alpha = 2;
			beta = 2;
			d1 = 1;
			d2 = 1;
			timevec = [0, 0.25, 0.5, 1, 2, 3];
			
		case 6
			
			% set 6:
			mu = 1;
			yh = 1;
			ye = 0.03;
			ks = 0.8;
			ki = 2;
			ke = 0;
			kl = 0.1;
			alpha = 2;
			beta = 2;
			d1 = 1;
			d2 = 0.0001;
			timevec = [0, 3, 6, 9, 12, 15];
			
	end
	
	% simplify parameters
	a1 = mu/yh;
	a2 = kl*ks;
	a3 = ke*ks;
	a4 = (kl+ki)*ks;
	a5 = ye*mu;
	
	% perform test to ensure strict diagonal dominance inside boundary
	while((dt*max([a2+a3, mu, a4, a5]) >= ks) && (dt*(mu + a5) >= ks))
		
		dt = dt/2;
		Rx = dt/dx/dx;
		Ry = dt/dy/dy;
		
	end
	
	% active biomass
	Li = 6;
	C = [0.25, 0.325, 0.275, 0.3, 0.2, 0.225];
	r = [100, 50, 30, 80, 90, 100];
	xy = [0.25, 0.3; 0.5, 0.25; 0.7, 0.65; 0.4, 0.8; 0.5, 0.55; ...
		0.8, 0.3];
	INITX = rinitprofile(XX, YY, Li,  C, r, xy)';
	
	% inert biomass
	Li = 5;
	C = [ 0.025, 0.03,  0.035, 0.2,  0.025];
	r = [25, 50, 125, 100, 50];
	xy = [0.3, 0.3; 0.55, 0.25; 0.76, 0.65; 0.45, 0.8; 0.55, 0.55];
	INITI = initprofile(XX, YY, Li,  C, r, xy)';
	
	% extracellular polymeric matrix
	Li = 4;
	C = [0.3, 0.35, 0.4, 0.25];
	r = [200, 300, 600, 400];
	xy = [0.2, 0.2; 0.8, 0.8; 0.3, 0.7; 0.7, 0.2];
	INITE = initprofile(XX, YY, Li,  C, r, xy)';
	
	% default initialization starting at t = 0
	if(p == 1 || p == 4)
		
		INITS = ones(size(XX))';
		
	elseif(p == 2 || p == 5)
		
		INITS = 0.2*ones(size(XX))';
		INITE = 0.4*ones(size(XX))';
		
	elseif(p == 3 || p == 6)
		
		INITS = exp(-5*(XX.^2+YY.^2))';
		INITX = INITX;
		INITI = zeros(size(XX))';
		INITE = zeros(size(XX))';
		
	end
	
	S = INITS';
	X = INITX';
	I = INITI';
	E = INITE';
	
	St = S';
	Xt = X';
	It = I';
	Et = E';
	
	
	% test if initial mass is bounded by 1
	if(X + I + E >= 1)
		
		break
		
	end
	
	% establish initial boundary conditions
	Bcond = ones(M+1, N+1);
	Bcond(:, 1) = 0;
	Bcond(:, N+1) = 0;
	Bcond(1, :) = 0;
	Bcond(M+1, :) = 0;
	
	if(Neu == 0)
		
		% for Dirichlet we fix substrate to walls at 1
		S(:, 1) = 1;
		S(:, M+1) = 1;
		S(1, :) = 1;
		S(N+1, 1) = 1;
		
		% for Dirichlet we fix biomass to ground at 0
		X = X.*Bcond';
		I = I.*Bcond';
		E = E.*Bcond';
		
	else
		
		% for Neumann we fix flow at boundaries to 0
		S = S.*Bcond';
		X = X.*Bcond';
		I = I.*Bcond';
		E = E.*Bcond';
		
	end
	
	% create a vector of all state variables
	vold = [reshape(S, (M+1)*(N+1), 1); reshape(X, (M+1)*(N+1), 1); ...
		reshape(I, (M+1)*(N+1), 1); reshape(E, (M+1)*(N+1), 1)];
	
	% preallocate a size(N+1, N+1) zero matrix
	Balloc = zeros(N+1);
	Z1 = Balloc;
	Z2 = Balloc;
	H = Balloc;
	G1 = Balloc;
	G2 = Balloc;
	BS = sparse(blkdiag(eye(N+1), zeros((M-1)*(N+1)), eye(N+1)));
	
	% associate Neuman or Dirichlet boundary conditions outside of loop
	for k = 1:N+1
		
		BS(k, N+1+k) = -Neu;
		BS((M+1)*(N+1)-k+1, (N+1)*(M+1)-k-N) = -Neu;
		
	end
	
	BS(1, 1) = 1+Neu;
	BS(1, 2) = -Neu;
	BS(N+1, N) = -Neu;
	BS(N+1, N+1) = 1+Neu;
	BS(M*(N+1) + 1, M*(N+1) + 1) = 1+Neu;
	BS(M*(N+1) + 1, M*(N+1) + 2) = -Neu;
	BS((M+1)*(N+1), M*(N+1) + N) = -Neu;
	BS((M+1)*(N+1), (M+1)*(N+1)) = 1+Neu;
	
	BX = BS;
	BI = BS;
	BE = BS;
	
	% commands for animation
	if animate == 1
		
		if(loop == 1)
			
			% set up the movie
			writerObj = VideoWriter('D.avi', 'Uncompressed AVI');
			writerObj.FrameRate = 1;
			myVideo.Quality = 100;
			open(writerObj);
			
		else
			
			writerObj = VideoWriter('N.avi', 'Uncompressed AVI');
			writerObj.FrameRate = 1;
			myVideo.Quality = 100;
			open(writerObj);
			
		end
		
	end
	
	if(waitbarflag == 1)
		
		h = waitbar(0,'Please wait...');
		steps = length(0:dt:T);
		step = 1;
		
	end
	
	for t = 0:dt:T
		
		% enforce solution to be bounded by 1
		%while(1 - dt*max(a5/ks, kl + max(ke, ki)) <= max(vold))
		%	
		%	dt = dt/2;
	%		Rx = dt/dx/dx;
	%		Ry = dt/dy/dy;
	%		
	%	end
		
		S = St;
		X = Xt;
		I = It;
		E = Et;
		
		% basic definitions for matrix assembly
		Alphapx = -Rx*Dfunc((X(3:M+1, 2:N) + X(2:M, 2:N) + ...
			I(3:M+1, 2:N) + I(2:M, 2:N) + E(3:M+1, 2:N) + ...
			E(2:M, 2:N))/2, alpha, beta, d2);
		Alphamx = -Rx*Dfunc((X(1:M-1, 2:N) + X(2:M, 2:N) + ...
			I(1:M-1, 2:N) + I(2:M, 2:N) + E(1:M-1, 2:N) + ...
			E(2:M, 2:N))/2, alpha, beta, d2);
		Alphapy = -Ry*Dfunc((X(2:M, 3:N+1) + X(2:M, 2:N) + ...
			I(2:M, 3:N+1) + I(2:M, 2:N) + E(2:M, 3:N+1) + ...
			E(2:M, 2:N))/2, alpha, beta, d2);
		Alphamy = -Ry*Dfunc((X(2:M, 1:N-1) + X(2:M, 2:N) + ...
			I(2:M, 1:N-1) + I(2:M, 2:N) + E(2:M, 1:N-1) + ...
			E(2:M, 2:N))/2, alpha, beta, d2);
		
		salpha = 1 - Alphamx - Alphapx - Alphamy - Alphapy;
		monod = dt./(ks + S(2:M,2:N));
		
		BetaS = salpha + a1*X(2:M, 2:N).*monod;
		BetaX = salpha + a4.*monod;
		BetaI = salpha;
		BetaE = salpha + a3.*monod;
		
		Zeta1 = -a2.*monod;
		Eta = -a3.*monod;
		Gamma1 = -mu*X(2:M,2:N).*monod;
		Zeta2 = -a4.*monod;
		Gamma2 = -a5*X(2:M,2:N).*monod;
		
		for m = 1:M-1
			
			% create block matrices
			BS((N+1)*m+1:(m+1)*(N+1), m*(N+1)-N:(m+2)*(N+1)) = ...
				[diag([0, Alphamx(m, :), 0]) ...
				eye(N+1)*tridiag([Alphamy(m, :), 0], ...
				[1, BetaS(m, :), 1], [0, Alphapy(m, :)]) ...
				diag([0, Alphapx(m, :), 0])];
			BX((N+1)*m+1:(m+1)*(N+1), m*(N+1)-N:(m+2)*(N+1)) = ...
				[diag([0, Alphamx(m, :), 0]) ...
				eye(N+1)*tridiag([Alphamy(m, :), 0], ...
				[1, BetaX(m, :), 1], [0, Alphapy(m, :)]) ...
				diag([0, Alphapx(m, :), 0])];
			BI((N+1)*m+1:(m+1)*(N+1), m*(N+1)-N:(m+2)*(N+1)) = ...
				[diag([0, Alphamx(m, :), 0]) ...
				eye(N+1)*tridiag([Alphamy(m, :), 0], ...
				[1, BetaI(m, :), 1], [0, Alphapy(m, :)]) ...
				diag([0, Alphapx(m, :), 0])];
			BE((N+1)*m+1:(m+1)*(N+1), m*(N+1)-N:(m+2)*(N+1)) = ...
				[diag([0, Alphamx(m, :), 0]) ...
				eye(N+1)*tridiag([Alphamy(m, :), 0], ...
				[1, BetaE(m, :), 1], [0, Alphapy(m, :)]) ...
				sparse(diag([0, Alphapx(m, :), 0]))];
			Z1 = sparse(blkdiag(Z1, diag([0, Zeta1(m, :), 0])));
			Z2 = sparse(blkdiag(Z2, diag([0, Zeta2(m, :), 0])));
			G1 = sparse(blkdiag(G1, diag([0, Gamma1(m, :), 0])));
			G2 = sparse(blkdiag(G2, diag([0, Gamma2(m, :), 0])));
			H = sparse(blkdiag(H, diag([0, Eta(m, :), 0])));
			
			% input boundary conditions
			BS(m*(N+1) + 1, m*(N+1) + 2) = -Neu;
			BX(m*(N+1) + 1, m*(N+1) + 2) = -Neu;
			BI(m*(N+1) + 1, m*(N+1) + 2) = -Neu;
			BE(m*(N+1) + 1, m*(N+1) + 2) = -Neu;
			BS((m+1)*(N+1), m*(N+1) + N) = -Neu;
			BX((m+1)*(N+1), m*(N+1) + N) = -Neu;
			BI((m+1)*(N+1), m*(N+1) + N) = -Neu;
			BE((m+1)*(N+1), m*(N+1) + N) = -Neu;
			
		end
		
		Z1 = sparse(blkdiag(Z1, zeros(N+1)));
		Z2 = sparse(blkdiag(Z2, zeros(N+1)));
		G1 = sparse(blkdiag(G1, zeros(N+1)));
		G2 = sparse(blkdiag(G2, zeros(N+1)));
		H = sparse(blkdiag(H, zeros(N+1)));
		
		% create M-matrix
		MM = sparse([BS Z1 zeros((M+1)*(N+1)) H;...
			G1 BX zeros((M+1)*(N+1)) zeros((M+1)*(N+1));...
			zeros((M+1)*(N+1)) Z2 BI zeros((M+1)*(N+1));...
			G2 zeros((M+1)*(N+1)) zeros((M+1)*(N+1)) BE]);
		
		% call solver to find next state vnew
		%bcgstest3;
		vnew = mldivide(MM, vold);
		lv = length(vnew);
		
		% revert back to matrices from vector form
		S = reshape(vnew(1:lv/4), N+1, M+1)';
		X = reshape(vnew(lv/4+1:lv/2), N+1, M+1)';
		I = reshape(vnew(lv/2+1:3*lv/4), N+1, M+1)';
		E = reshape(vnew(3*lv/4+1:lv), N+1, M+1)';
		
		% save results
		St = S;
		Xt = X;
		It = I;
		Et = E;
		
		% reset boundary conditions
		if(Neu == 0)
			
			% for Dirichlet we fix substrate to walls
			S(:, 1) = 1;
			S(:, N+1) = 1;
			S(1, :) = 1;
			S(M+1, 1) = 1;
			
			% for Dirichlet we fix biomass to ground
			X = X.*Bcond;
			I = I.*Bcond;
			E = E.*Bcond;
			
		else
			
			% for Neumann we fix zero flow at boundaries
			S = S.*Bcond;
			X = X.*Bcond;
			I = I.*Bcond;
			E = E.*Bcond;
			
		end
		
		vold = [reshape(S', (M+1)*(N+1), 1); ...
			reshape(X', (M+1)*(N+1), 1); ...
			reshape(I', (M+1)*(N+1), 1); reshape(E', (M+1)*(N+1), 1)];
		
		Z1 = Balloc;
		Z2 = Balloc;
		H = Balloc;
		G1 = Balloc;
		G2 = Balloc;
		BS = sparse(blkdiag(eye(N+1), zeros((M-1)*(N+1)), eye(N+1)));
		
		% associate Neuman or Dirichlet boundary conditions
		for k = 1:N+1
			
			BS(k, N+1+k) = -Neu;
			BS((M+1)*(N+1)-k+1, (N+1)*(M+1)-k-N) = -Neu;
			
		end
		
		BS(1, 1) = 1+Neu;
		BS(1, 2) = -Neu;
		BS(N+1, N) = -Neu;
		BS(N+1, N+1) = 1+Neu;
		BS(M*(N+1) + 1, M*(N+1) + 1) = 1+Neu;
		BS(M*(N+1) + 1, M*(N+1) + 2) = -Neu;
		BS((M+1)*(N+1), M*(N+1) + N) = -Neu;
		BS((M+1)*(N+1), (M+1)*(N+1)) = 1+Neu;
		
		BX = BS;
		BI = BS;
		BE = BS;
		
		if(animate == 1)
			
			surfinset3;
			frame = getframe(gcf);
			writeVideo(writerObj, frame);
			% pause(2);
			closereq;
			
		end
		
		% for current memory used flag this
		if checkmemory == 1
			
			memory;
			
		end
		
		if(waitbarflag == 1)
			
			step = step + 1;
			waitbar(step / steps, h);
			
		end
	
		if(animate ~= 1)
			
			if(Neu == 0 &&  any(abs(t - timevec) < 1e-8))
				
				dird = 'C:\Users\Richard\Desktop\Dissertation\DE3';
				file = [dird, '_', num2str(t), '_', num2str(p)];
				save([file, '.mat']);
				surfinset2;
				export_fig(file,'-png','-jpg','-tiff', '-transparent');
				
			elseif(Neu == 1 && any(abs(t - timevec) < 1e-8))
				
				dirn = 'C:\Users\Richard\Desktop\Dissertation\NE3';
				file = [dirn, '_', num2str(t), '_', num2str(p)];
				save([file, '.mat']);
				%surfinset2;
				export_fig(file,'-png','-jpg','-tiff', '-transparent');
				
			end
						
		else
			
			close(writerObj);
			
		end
		
	end
	
	if(waitbarflag == 1)
		
		close(h);
		
	end
	
	testresidual;
	
	% give total execution time for entire test run
	toc;	
	
	if(loop == 1)
		
		UX = Xt;
		US = St;
		UI = It;
		UE = Et;
		
	elseif(loop == 2)
		
		VX = Xt;
		VS = St;
		VI = It;
		VE = Et;
		
	end
	
end

if(animate ~= 1)
	
	comparecontour3(XX, YY, UX, VX, US, VS, UI, VI, UE, VE)
	
end