% once dirichletbio2d runs to completion, max values
% are set and this program standardizes the axes

keep animate paramcond tmaxS tmaxX tmaxE tmaxI INITS INITX INITI INITE dx dy dt T
closereq
tic;


% establish parameters
Rx = dt/dx/dx;
Ry = dt/dy/dy;

switch paramcond
    case 1
        %SET 1;
        mu = 0.002;
        yh = 0.03;
        kl = 0.1;
        ks = 0.2;
        ke = 0.03;
        ki = 0.04;
        ye = 0.03;
        alpha = 2;
        beta = 2;
        %clear maxS maxX maxI maxE;
        
    case 2
        %SET 2:
        mu = 0.4;
        yh = 0.03;
        ye = 0.03;
        ks = 0.2;
        ki = 0.0004;
        ke = 0.03;
        kl = 0.0001;
        alpha = 4;
        beta  = 4;
        %clear maxS maxX maxI maxE;
        
    case 3
        %SET 3:
        mu = 2;
        yh = 0.35;
        ye = 0.03;
        ks = 0.2;
        ki = 0.4;
        ke = 0.03;
        kl = 0.0001;
        alpha = 4;
        beta = 4;
        %clear maxS maxX maxI maxE;
        
end

a1 = mu/yh;
a2 = kl*ks;
a3 = ke*ks;
a4 = (kl+ki)*ks;
a5 = ye*mu;

x = 0:dx:1;
y = 0:dy:1;
M = length(x)-1;
N = length(y)-1;
[XX, YY] = meshgrid(x, y);

S = INITS';
X = INITX';
I = INITI';
E = INITE';

% establish initial Dirichlet boundary conditions
S(1:M+1, 1) = 0;
S(1:M+1, N+1) = 0;
X(1:M+1, 1) = 0;
X(1:M+1, N+1) = 0;
I(1:M+1, 1) = 0;
I(1:M+1, N+1) = 0;
E(1:M+1, 1) = 0;
E(1:M+1, N+1) = 0;
S(1, 1:N+1) = 0;
S(M+1, 1:N+1) = 0;
X(1, 1:N+1) = 0;
X(M+1, 1:N+1) = 0;
I(1, 1:N+1) = 0;
I(M+1, 1:N+1) = 0;
E(1, 1:N+1) = 0;
E(M+1, 1:N+1) = 0;

% create a vector of all state variables
b0 = [reshape(S', (M+1)*(N+1), 1); reshape(X', (M+1)*(N+1), 1); reshape(I', (M+1)*(N+1), 1); reshape(E', (M+1)*(N+1), 1)];
b = b0;
Z1 = zeros(N+1);
Z2 = zeros(N+1);
H = zeros(N+1);
G1 = zeros(N+1);
G2 = zeros(N+1);

BS = sparse(blkdiag(eye(N+1), zeros((M-1)*(N+1)), eye(N+1)));
BX = sparse(blkdiag(eye(N+1), zeros((M-1)*(N+1)), eye(N+1)));
BI = sparse(blkdiag(eye(N+1), zeros((M-1)*(N+1)), eye(N+1)));
BE = sparse(blkdiag(eye(N+1), zeros((M-1)*(N+1)), eye(N+1)));

movieiterate = 1;
for t=0:dt:T
    
    Alphapx = -Rx*D((X(3:M+1, 2:N) + X(2:M, 2:N) + I(3:M+1, 2:N) + I(2:M, 2:N) + E(3:M+1, 2:N) + E(2:M, 2:N))/2, alpha, beta);
    Alphamx = -Rx*D((X(1:M-1, 2:N) + X(2:M, 2:N) + I(1:M-1, 2:N) + I(2:M, 2:N) + E(1:M-1, 2:N) + E(2:M, 2:N))/2, alpha, beta);
    Alphapy = -Ry*D((X(2:M, 3:N+1) + X(2:M, 2:N) + I(2:M, 3:N+1) + I(2:M, 2:N) + E(2:M, 3:N+1) + E(2:M, 2:N))/2, alpha, beta);
    Alphamy = -Ry*D((X(2:M, 1:N-1) + X(2:M, 2:N) + I(2:M, 1:N-1) + I(2:M, 2:N) + E(2:M, 1:N-1) + E(2:M, 2:N))/2, alpha, beta);
    BetaS = 1 - Alphamx - Alphapx - Alphamy - Alphapy + a1*dt*X(2:M, 2:N)./(ks + S(2:M,2:N));
    BetaX = 1 - Alphamx - Alphapx - Alphamy - Alphapy + a4*dt./(ks + S(2:M,2:N));
    BetaI = 1 - Alphamx - Alphapx - Alphamy - Alphapy;
    BetaE = 1 - Alphamx - Alphapx - Alphamy - Alphapy + a3*dt./(ks + S(2:M,2:N));
    
    Zeta1 = -a2*dt./(ks + S(2:M,2:N));
    Eta = -a3*dt./(ks + S(2:M,2:N));
    Gamma1 = -mu*dt*X(2:M,2:N)./(ks + S(2:M,2:N));
    Zeta2 = -a4*dt./(ks + S(2:M,2:N));
    Gamma2 = -a5*dt*X(2:M,2:N)./(ks + S(2:M,2:N));
    
    for m = 1:M-1
        BS((N+1)*m+1:(m+1)*(N+1), m*(N+1)-N:(m+2)*(N+1)) = [sparse(diag([0, Alphamx(m, :), 0])) eye(N+1)*tridiag([Alphamy(m, :), 0], [1, BetaS(m, :), 1], [0, Alphapy(m, :)]) sparse(diag([0, Alphapx(m, :), 0]))];
        BX((N+1)*m+1:(m+1)*(N+1), m*(N+1)-N:(m+2)*(N+1)) = [sparse(diag([0, Alphamx(m, :), 0])) eye(N+1)*tridiag([Alphamy(m, :), 0], [1, BetaX(m, :), 1], [0, Alphapy(m, :)]) sparse(diag([0, Alphapx(m, :), 0]))];
        BI((N+1)*m+1:(m+1)*(N+1), m*(N+1)-N:(m+2)*(N+1)) = [sparse(diag([0, Alphamx(m, :), 0])) eye(N+1)*tridiag([Alphamy(m, :), 0], [1, BetaI(m, :), 1], [0, Alphapy(m, :)]) sparse(diag([0, Alphapx(m, :), 0]))];
        BE((N+1)*m+1:(m+1)*(N+1), m*(N+1)-N:(m+2)*(N+1)) = [sparse(diag([0, Alphamx(m, :), 0])) eye(N+1)*tridiag([Alphamy(m, :), 0], [1, BetaE(m, :), 1], [0, Alphapy(m, :)]) sparse(diag([0, Alphapx(m, :), 0]))];
        Z1 = sparse(blkdiag(Z1, diag([0, Zeta1(m, :), 0])));
        Z2 = sparse(blkdiag(Z2, diag([0, Zeta2(m, :), 0])));
        G1 = sparse(blkdiag(G1, diag([0, Gamma1(m, :), 0])));
        G2 = sparse(blkdiag(G2, diag([0, Gamma2(m, :), 0])));
        H = sparse(blkdiag(H, diag([0, Eta(m, :), 0])));
    end
    
    Z1 = sparse(blkdiag(Z1, zeros(N+1)));
    Z2 = sparse(blkdiag(Z2, zeros(N+1)));
    G1 = sparse(blkdiag(G1, zeros(N+1)));
    G2 = sparse(blkdiag(G2, zeros(N+1)));
    
    H = blkdiag(H, zeros(N+1));
    Hsparse = sparse(H);
    clear H;
    MM = [BS Z1 zeros((M+1)*(N+1)) Hsparse;...
        G1 BX zeros((M+1)*(N+1)) zeros((M+1)*(N+1));...
        zeros((M+1)*(N+1)) Z2 BI zeros((M+1)*(N+1));...
        G2 zeros((M+1)*(N+1)) zeros((M+1)*(N+1)) BE];
    MMsparse = sparse(MM);
    clear MM;
    b = bicgstab(MMsparse, b, 1e-6, 200);
    
    S = reshape(b(1:length(b)/4), N+1, M+1)';
    X = reshape(b(length(b)/4+1:2*length(b)/4), N+1, M+1)';
    I = reshape(b(2*length(b)/4+1:3*length(b)/4), N+1, M+1)';
    E = reshape(b(3*length(b)/4+1:length(b)), N+1, M+1)';
    
    S(1:M+1, 1) = 0;
    S(1:M+1, N+1) = 0;
    X(1:M+1, 1) = 0;
    X(1:M+1, N+1) = 0;
    I(1:M+1, 1) = 0;
    I(1:M+1, N+1) = 0;
    E(1:M+1, 1) = 0;
    E(1:M+1, N+1) = 0;
    S(1, 1:N+1) = 0;
    S(M+1, 1:N+1) = 0;
    X(1, 1:N+1) = 0;
    X(M+1, 1:N+1) = 0;
    I(1, 1:N+1) = 0;
    I(M+1, 1:N+1) = 0;
    E(1, 1:N+1) = 0;
    E(M+1, 1:N+1) = 0;
    Z1 = zeros((N+1));
    Z2 = zeros((N+1));
    H = zeros((N+1));
    G1 = zeros((N+1));
    G2 = zeros((N+1));
    
    
    if(t == 0 || t == 0.1 || t == 0.2 || t == 0.5 || t == 1 || t == 2 || t == 5)
        %surfacegraph
        insetplot
    end
    
    if animate == 1
        %maximize
        fprintf('Frame: %03d\n', t);
        moviefilename = strcat('animation', num2str(paramcond), '_', num2str(T), '.avi');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % insert plot here
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        hf = figure;
        grid off;
        set(hf,'position',[10 190 640 500]);
        subplot(2, 2, 1)
        surf(XX,YY,S','FaceColor','interp',...
            'EdgeColor','k',...
            'FaceLighting','phong');
        title(strcat('Surface Plot of S at T = ', num2str(t)));%,'Interpreter','LaTex');
        xlabel('x');%,'Interpreter','LaTex');
        ylabel('y');%,'Interpreter','LaTex');
        zlabel(strcat('S(x, y, ', num2str(t), ')'));%,'Interpreter','LaTex');
        %daspect([5 5 0.1]);
        axis([0, 1, 0, 1, 0, tmaxS]);
        view(-50,30);
        camlight left;
        
        
        grid off;
        subplot(2, 2, 2)
        surf(XX,YY,X','FaceColor','interp',...
            'EdgeColor','k',...
            'FaceLighting','phong');
        title(strcat('Surface Plot of X at T = ', num2str(t)));%,'Interpreter','LaTex');
        xlabel('x');%,'Interpreter','LaTex');
        ylabel('y');%,'Interpreter','LaTex');
        zlabel(strcat('X(x, y, ', num2str(t), ')'));%,'Interpreter','LaTex');
        %daspect([5 5 0.1]);
        axis([0,1, 0, 1, 0, tmaxX]);
        view(-50,30);
        camlight left;
        
        
        grid off;
        subplot(2, 2, 3)
        surf(XX,YY,I','FaceColor','interp',...
            'EdgeColor','k',...
            'FaceLighting','phong');
        title(strcat('Surface Plot of I at T = ', num2str(t)));%,'Interpreter','LaTex');
        xlabel('x');%,'Interpreter','LaTex');
        ylabel('y');%,'Interpreter','LaTex');
        zlabel(strcat('I(x, y, ', num2str(t), ')'));%,'Interpreter','LaTex');
        %daspect([5 5 0.1]);
        axis([0,1, 0, 1, 0, tmaxI]);
        view(-50,30);
        camlight left;
        
        grid off;
        subplot(2, 2, 4)
        surf(XX,YY,E','FaceColor','interp',...
            'EdgeColor','k',...
            'FaceLighting','phong');colormap(hot)
        title(strcat('Surface Plot of E at T = ', num2str(t)));%,'Interpreter','LaTex');
        xlabel('x');%,'Interpreter','LaTex');
        ylabel('y');%,'Interpreter','LaTex');
        zlabel(strcat('E(x, y, ', num2str(t), ')'));%,'Interpreter','LaTex');
        %daspect([5 5 0.1]);
        axis([0,1, 0, 1, 0, tmaxE]);
        view(-50,30);
        camlight left;
        grid off
        mov1(movieiterate) = getframe(hf);
        pause(0.05);
        closereq
        movieiterate = movieiterate + 1;
    end
    
end

if animate == 1
    movie2avi(mov1, moviefilename);
    set(gcf,'position',[10 190 640 500]);
    movie(gcf, mov1);
end

