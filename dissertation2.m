%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Biofilm Experiment 2                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning off; %#ok<*WNOFF>
clc;
clear all;

% begin recording time of total program
tic;

% end time
T = 1.01;

% evaluate using predefined paramaters and initial conditions
p = 2;

% display total amount of memory used
checkmemory = 0;

% flag to create and save animation and not save images
animate = 1;

% flag to create a window showing estimated time left
waitbarflag = 0;

% loop through two conditions only changing one setting
for loop = 2:2
	u = 0;
	
	% establish intial spatial and temporal resolution
	dx = 0.01;
	dt = 0.01;
	
	% this performs loop over Neumann and Dirichlet conditions
	if(loop == 1)
		
		Neu = 0;
		
	else
		
		Neu = 1;
		
	end
	
	% mesh grid settings, number of node points
	Rx = dt/dx/dx;
	x = 0:dx:1;
	t = 0:dt:T;
	M = length(x)-1;
	
	% predefined settings for parameters
	switch p
		
		case 1
			
			% set 1:
			mu = 1;
			yh = 1;
			kl = 0.01;
			ks = 0.008;
			ke = 0.03;
			ki = 0.4;
			ye = 0.03;
			alpha = 4;
			beta = 4;
			d1 = 1;
			d2 = 1;
			timevec = [0, 0.25, 0.5, 0.75, 1, 1.5, 2, 5, 10];
			
		case 2
			
			% set 2:
			mu = 1;
			yh = 1;
			kl = 0.03;
			ks = 0.8;
			ke = 0.03;
			ki = 0.04;
			ye = 0.03;
			alpha = 6;
			beta = 1;
			d1 = 1;
			d2 = 0.001;
			timevec = [0, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, ...
				9, 10, 12, 14, 16, 18, 20, 22, 24. 26. 28 30, 35, 40, ...
				45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 103];
			
		case 3
			
			% set 3:
			mu = 2;
			yh = 5;%0.35
			ye = 0.03;
			ks = .99;
			ki = 0.4;
			ke = 0.03;
			kl = 0.0001;
			alpha = 4;
			beta = 4;
			d1 = 1;
			d2 = 1;
			timevec = [0, 0.25, 0.5, 0.75, 1, 1.5, 2, 5, 10];
			
	end
	
	% simplify parameters
	a1 = mu/yh;
	a2 = kl*ks;
	a3 = ke*ks;
	a4 = (kl+ki)*ks;
	a5 = ye*mu;
	
	% perform test to ensure strict diagonal dominance inside boundary
	%while((dt*max([a2+a3, mu, a4, a5]) >= ks) || (dt*(mu + a5) >= ks))
	%
	%	dt = dt/2;
	%		Rx = dt/dx/dx;
	
	%	end
	
	% substrate
	Li = 1;
	C = [0.8];
	r = [10];
	x0 = [0.0];
	INITS = initprofile2(x, Li, C, r, x0);
	
	% active biomass
	Li = 2;
	C = [0.15, 0.2];
	r = [200, 150];
	x0 = [0.3, 0.65];
	INITX = initprofile2(x, Li, C, r, x0);
	
	% inert biomass
	Li = 2;
	C = [0.02,  0.01];
	r = [80, 60];
	x0 = [0.0, 1.0];
	INITI = initprofile2(x, Li, C, r, x0);
	
	% extracellular polymeric matrix
	Li = 6;
	C = 0.5*[0.1, 0.3, 0.35, 0.38, 0.4, 0.25];
	r = [100, 200, 300, 500, 600, 400];
	x0 = [0.0, 0.2, 0.45, 0.8, 0.3, 0.7];
	INITE = initprofile2(x, Li, C, r, x0);
	%INITE = 0.05*ones(size(x))
	
	% default initialization starting at t = 0
	S = INITS';
	X = INITX';
	I = INITI';
	E = INITE';
	
	St = S;
	Xt = X;
	It = I;
	Et = E;
	
	% test if initial mass is bounded by 1
	if(X + I + E >= 1)
		
		break
		
	end
	
	if(Neu == 0)
		
		% for Dirichlet we fix substrate to walls at initial
		S0 = S(1);
		SL = S(M+1);
		
		% for Dirichlet we fix biomass to ground at 0
		X(1) = 0;
		X(M+1) = 0;
		I(1) = 0;
		I(M+1) = 0;
		E(1) = 0;
		E(M+1) = 0;
		
	else
		
		% for Neumann we fix flow at boundaries to 0
		S(1) = 0;
		S(M+1) = 0;
		X(1) = 0;
		X(M+1) = 0;
		I(1) = 0;
		I(M+1) = 0;
		E(1) = 0;
		E(M+1) = 0;
		
	end
	
	% create a vector of all state variables
	vold = [S' X' I' E'];
	
	% preallocate a size(M+1, M+1) zero matrix
	Balloc = sparse(zeros(M+1));
	Z1 = Balloc;
	Z2 = Balloc;
	H = Balloc;
	G1 = Balloc;
	G2 = Balloc;
	
	% commands for animation
	if(animate == 1)
		
		% set up the movie
		if(loop == 1)
			
			writerObj = VideoWriter('out1.avi', 'Uncompressed AVI');
			writerObj.FrameRate = 1; % How many frames per second.
			myVideo.Quality = 100;
			open(writerObj);
			
		else
			
			writerObj = VideoWriter('out2.avi', 'Uncompressed AVI');
			writerObj.FrameRate = 0.1; % How many frames per second.
			myVideo.Quality = 10;
			open(writerObj)
			
		end
		
	end
	
	if(waitbarflag == 1)
		
		h = waitbar(0,'Please wait...');
		steps = length(0:dt:T);
		step = 1;
		
	end
	
	SS = [];
	XX = [];
	II = [];
	EE = [];
	
	for t = 0:dt:T
		
		S = St;
		X = Xt;
		I = It;
		E = Et;
		
		if(mod(u, 10) == 0)
			t
			SS = [SS St];
			XX = [XX Xt];
			II = [II It];
			EE = [EE Et];
		end
		
		%% enforce solution to be bounded by 1
		%while(1 - dt*max(a5/ks, kl + max(ke, ki)) <= vold)
		
		%		dt = dt/2;
		%		Rx = dt/dx/dx;
		
		%	end
		
		% basic definitions for matrix assembly
		Alphapx = -Rx*Dfunc((X(3:M+1) + X(2:M) + ...
			I(3:M+1) + I(2:M) + E(3:M+1) + ...
			E(2:M))/2, alpha, beta, d2);
		Alphamx = -Rx*Dfunc((X(1:M-1) + X(2:M) + ...
			I(1:M-1) + I(2:M) + E(1:M-1) + ...
			E(2:M))/2, alpha, beta, d2);
		
		salpha = 1 - Alphamx - Alphapx;
		monod = dt./(ks + S(2:M));
		
		BetaS = salpha + a1*X(2:M).*monod;
		BetaX = salpha + a4.*monod;
		BetaI = salpha;
		BetaE = salpha + a3.*monod;
		
		Zeta1 = -a2.*monod;
		Zeta2 = -a4.*monod;
		
		Eta = -a3.*monod;
		Gamma1 = -mu*X(2:M).*monod;
		Gamma2 = -a5*X(2:M).*monod;
		
		% create block matrices and input boundary conditions
		BS = tridiag([Alphamx; -Neu], [1; BetaS; 1], [-Neu; Alphapx]);
		BX = tridiag([Alphamx; -Neu], [1; BetaX; 1], [-Neu; Alphapx]);
		BI = tridiag([Alphamx; -Neu], [1; BetaI; 1], [-Neu; Alphapx]);
		BE = tridiag([Alphamx; -Neu], [1; BetaE; 1], [-Neu; Alphapx]);
		
		Z1 = sparse(diag([0; Zeta1; 0]));
		Z2 = sparse(diag([0; Zeta2; 0]));
		G1 = sparse(diag([0; Gamma1; 0]));
		G2 = sparse(diag([0; Gamma2; 0]));
		H = sparse(diag([0; Eta; 0]));
		
		% create M-matrix
		MM = sparse([BS Z1 Balloc H;...
			G1 BX Balloc Balloc;...
			Balloc Z2 BI Balloc;...
			G2 Balloc Balloc BE]);
		
		% call solver to find next state vnew
		bcgstest2;
		lv = length(vnew);
		
		% revert back to individual state variables from vector form
		S = vnew(1:lv/4);
		X = vnew(lv/4+1:lv/2);
		I = vnew(lv/2+1:3*lv/4);
		E = vnew(3*lv/4+1:lv);
		
		% save results for LHS
		St = S;
		Xt = X;
		It = I;
		Et = E;
		
		if(Neu == 0)
			
			% for Dirichlet we fix substrate to walls at initial
			S(1) = S0;
			S(M+1) = SL;
			
			% for Dirichlet we fix biomass to ground at 0
			X(1) = 0;
			X(M+1) = 0;
			I(1) = 0;
			I(M+1) = 0;
			E(1) = 0;
			E(M+1) = 0;
			
		else
			
			% for Neumann we fix flow at boundaries to 0
			S(1) = 0;
			S(M+1) = 0;
			X(1) = 0;
			X(M+1) = 0;
			I(1) = 0;
			I(M+1) = 0;
			E(1) = 0;
			E(M+1) = 0;
			
		end
		
		vold = [reshape(S', M+1, 1); reshape(X', M+1, 1); ...
			reshape(I', M+1, 1); reshape(E', M+1, 1)]';
		
		Z1 = Balloc;
		Z2 = Balloc;
		H = Balloc;
		G1 = Balloc;
		G2 = Balloc;
		if(animate == 1)
			
			surfinset2;
			frame = getframe(gcf);
			writeVideo(writerObj, frame);
			%pause(2);
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
		
		% plot commands if animation is not used
		if(animate ~= 1)
			
			if(Neu == 0 &&  any(abs(t - timevec) < 1e-8))
				
				dird = 'C:\Users\Richard\Desktop\Dissertation\DE2';
				file = [dird, '_', num2str(t), '_', num2str(p)];
				save([file, '.mat']);
				surfinset2;
				export_fig(file,'-png','-jpg','-tiff', '-transparent');
				
			elseif(Neu == 1 && any(abs(t - timevec) < 1e-8))
				
				dirn = 'C:\Users\Richard\Desktop\Dissertation\NE2';
				file = [dirn, '_', num2str(t), '_', num2str(p)];
				save([file, '.mat']);
				surfinset2;
				export_fig(file,'-png','-jpg','-tiff', '-transparent');
				
			end
			
		else
			
			%close(writerObj);
			
		end
		
		u = u + 1;
	end
	
	if(animate == 1)
		
		close(writerObj);
		
	end
	
	if(waitbarflag == 1)
		
		close(h);
		
	end
	
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
	
	figure;
	comparecontour2(x, UX, VX, US, VS, UI, VI, UE, VE)
	
end


