function u=new_approximation(b,k,M,N,Rx2,Ry2,A,K1,K2,K3,K4,Neu)
  A=matrix(M,N,Rx2,Ry2,k,b,A,K1,K2,K3,K4);
  for n=1:N
    b(n)=1-Neu;     b(M*N+1-n)=1-Neu; 
    b(M*N+n)=0; b(2*M*N+1-n)=0;
  end
  for m=1:M-2
    b(m*N+1)=1-Neu;     b((m+1)*N)=1-Neu;
    b((M+m)*N+1)=0; b((M+m+1)*N)=0;
  end
  u=bicgstab(A,b,1e-6,40);
return