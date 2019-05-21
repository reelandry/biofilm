%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dissertation1.m                                                         %
% Biofilm Experiment 1                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% start from previous run
prevflag = 0;

if(prevflag == 0)
	
	warning off; %#ok<*WNOFF>
	clc;
	keep prevflag;
	
end

% begin recording time of total program
tic;

% end time
T = 2.01;

% evaluate using predefined paramaters and initial conditions
p = 3;

% display total amount of memory used
checkmemory = 0;

% flag to create and save animation and not save images
animate = 0;

% flag to create a window showing estimated time left
waitbarflag = 1;

% loop through two conditions only changing one setting
for loop = 2:2
	
	% establish intial spatial and temporal resolution
	dx = 0.01;
	dy = 0.01;
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
	M = length(x)-1;
	N = length(y)-1;
	[XX, YY] = meshgrid(x, y);
	
	% predefined settings for parameters
	switch p
		
		case 1
			
			% set 1
			K1 = 0;
			K2 = 0;
			K3 = 0.4;
			K4 = 0;
			alpha = 4;
			beta = 4;
			d1 = 0;
			d2 = 0.0001;
			timevec = [0, 1, 2, 3, 4, 5, 6, 7, 8];
			
		case 2
			
			% set 2:
			K1 = 0.85;
			K2 = 0.0012;
			K3 = 0.4;
			K4 = 0.3;
			alpha = 4;
			beta = 4;
			d1 = 0.002;
			d2 = 0.0001;
			timevec = [0, 2, 4, 6, 8, 10, 12, 14];
			
		case 3
			
			% set 3:
			K1 = 0.65;
			K2 = 0.36;
			K3 = 0.2;
			K4 = 0.3;
			alpha = 4;
			beta = 4;
			d1 = 0.0015;
			d2 = 0.0001;
			timevec = [0, 1, 2.5, 5, 7.5 10, 15, 20, 30];
			
	end
	
	% perform test to ensure strict diagonal dominance inside boundary
	while(dt*(K3-K2) >= 1)
		
		dt = dt/2
		Rx = dt/dx/dx;
		Ry = dt/dy/dy;
		
	end
	
	if(prevflag == 0)
	% active biomass initial conditions
	Li = 5;
	C = [0.025, 0.03, 0.035, 0.02, 0.025];
	r = [25, 50, 125, 100, 50];
	xy = [0.25, 0.3; 0.5, 0.25; 0.7, 0.65; 0.4, 0.8; 0.5, 0.55];
	INITX = initprofile(XX, YY, Li, C, r, xy)';
	
	% default initialization starting at t = 0
	INITS = ones(size(XX))';
	S = INITS';
	X = INITX';
	
	St = S';
	Xt = X';
	
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
		
	else
		
		% for Neumann we fix flow at boundaries to 0
		S = S.*Bcond';
		X = X.*Bcond';
		
	end
	
	% create a vector of all state variables
	vold = [reshape(S, (M+1)*(N+1), 1); reshape(X, (M+1)*(N+1), 1)];
	lv = length(vold);
	
	% preallocate a size(N+1, N+1) zero matrix
	Balloc = zeros(N+1);
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
	
	end
	
	t = 14.02;
	
	% main loop for iterations over time
	while(t <= T)
		
		S = St;
		X = Xt;
					
		% basic definitions for matrix assembly
		Phipx = -Rx*d1*X(2:M, 2:N);
		Phimx = -Rx*d1*X(2:M, 2:N);
		Phipy = -Ry*d1*X(2:M, 2:N);
		Phimy = -Ry*d1*X(2:M, 2:N);
		
		Psipx = -Rx*Dfunc((X(3:M+1, 2:N) + X(2:M, 2:N)/2), ...
			alpha, beta, d2);
		Psimx = -Rx*Dfunc((X(1:M-1, 2:N) + X(2:M, 2:N)/2), ...
			alpha, beta, d2);
		Psipy = -Ry*Dfunc((X(2:M, 3:N+1) + X(2:M, 2:N)/2), ...
			alpha, beta, d2);
		Psimy = -Ry*Dfunc((X(2:M, 1:N-1) + X(2:M, 2:N)/2), ...
			alpha, beta, d2);
		
		sPhi = 1 - Phimx - Phipx - Phimy - Phipy;
		sPsi = 1 - Psimx - Psipx - Psimy - Psipy;
		monod = dt./(K4 + S(2:M,2:N));
		
		PhiS = sPhi + K1*X(2:M, 2:N).*monod;
		ChiX = sPsi + dt*K2 - K3*S(2:M, 2:N).*monod;
		
		for m = 1:M-1
			
			% create block matrices
			BS((N+1)*m+1:(m+1)*(N+1), m*(N+1)-N:(m+2)*(N+1)) = ...
				[diag([0, Phimx(m, :), 0]) ...
				eye(N+1)*tridiag([Phimy(m, :), 0], ...
				[1, PhiS(m, :), 1], [0, Phipy(m, :)]) ...
				diag([0, Phipx(m, :), 0])];
			BX((N+1)*m+1:(m+1)*(N+1), m*(N+1)-N:(m+2)*(N+1)) = ...
				[diag([0, Psimx(m, :), 0]) ...
				eye(N+1)*tridiag([Psimy(m, :), 0], ...
				[1, ChiX(m, :), 1], [0, Psipy(m, :)]) ...
				diag([0, Psipx(m, :), 0])];
			
			% input boundary conditions
			BS(m*(N+1) + 1, m*(N+1) + 2) = -Neu;
			BX(m*(N+1) + 1, m*(N+1) + 2) = -Neu;
			BS((m+1)*(N+1), m*(N+1) + N) = -Neu;
			BX((m+1)*(N+1), m*(N+1) + N) = -Neu;
			
		end
		
		% create M-matrix
		MM = sparse([BS zeros((M+1)*(N+1)); zeros((M+1)*(N+1)) BX]);
		
		% call solver to find next state vnew
		%bcgstest1;
		vworking = vold;
		vold = mldivide(MM, vold);
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% catch if solution is bad
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		% enforce solution to be bounded by 1
		voldu = vold(lv/2+1:lv);		
		
		
		% revert back to matrices from vector form
		S = reshape(vold(1:lv/2), N+1, M+1)';
		X = reshape(vold(lv/2+1:lv), N+1, M+1)';
		
		% save results
		St = S;
		Xt = X;
		
		% plot commands if animation is not used
		if(animate ~= 1)
			
			if(Neu == 0 &&  any(abs(t - timevec) < 1e-8))
				
				dird = 'C:\Users\Richard\Desktop\Dissertation\DE1';
				file = [dird, '_', num2str(t), '_', num2str(p)];
				save([file, '.mat']);
				surfinset1;
				export_fig(file,'-png','-jpg','-tiff', '-transparent');
				
			elseif(Neu == 1 && any(abs(t - timevec) < 1e-8))
				
				dirn = 'C:\Users\Richard\Desktop\Dissertation\NE1';
				file = [dirn, '_', num2str(t), '_', num2str(p)];
				save([file, '.mat']);
				surfinset1;
				export_fig(file,'-png','-jpg','-tiff', '-transparent');
				
			end
			
		else
			
			close(writerObj);
			
		end
		
		% reset boundary conditions
		if(Neu == 0)
			
			% for Dirichlet we fix substrate to walls
			S(:, 1) = 1;
			S(:, N+1) = 1;
			S(1, :) = 1;
			S(M+1, 1) = 1;
			
			% for Dirichlet we fix biomass to ground
			X = X.*Bcond;
			
		else
			
			% for Neumann we fix zero flow at boundaries
			S = S.*Bcond;
			X = X.*Bcond;
			
		end
		
		vold = [reshape(S', (M+1)*(N+1), 1); reshape(X', (M+1)*(N+1), 1)];
		
		BS = sparse(blkdiag(eye(N+1), zeros((M-1)*(N+1)), eye(N+1)));
		
		% reassociate Neumann or Dirichlet boundary conditions
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
		
		if(animate == 1)
			
			surfinset1;
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
		
		t = t + dt;
		
		
		
	end
	
	if(waitbarflag == 1)
		
		close(h);
		
	end
	
	% give total execution time for entire test run
	toc;
	
	if(loop == 1)
		
		US = St;
		UX = Xt;
		
	elseif(loop == 2)
		
		VS = St;
		VX = Xt;
		
	end
	
end

if(animate ~= 1)
	
	% compare the looped conditions
	comparecontour1(XX, YY, US, VS, UX, VX, M, N);
	
end

