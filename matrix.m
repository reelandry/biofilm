function A=matrix(M,N,Rx2,Ry2,k,u,A,K1,K2,K3,K4)
  for m=1:M-2
    for n=2:N-1
      A(m*N+n,m*N+n)=A(m*N+n,m*N+n)+K1*k*u((M+m)*N+n)/(K4+u(m*N+n));
      A((M+m)*N+n,(M+m)*N+n-1)=-Ry2*D(0.5*(u((M+m)*N+n-1)+u((M+m)*N+n)));
      A((M+m)*N+n,(M+m)*N+n+1)=-Ry2*D(0.5*(u((M+m)*N+n+1)+u((M+m)*N+n)));
      A((M+m)*N+n,(M+m-1)*N+n)=-Rx2*D(0.5*(u((M+m-1)*N+n)+u((M+m)*N+n)));
      A((M+m)*N+n,(M+m+1)*N+n)=-Rx2*D(0.5*(u((M+m+1)*N+n)+u((M+m)*N+n)));
      A((M+m)*N+n,(M+m)*N+n)=1-A((M+m)*N+n,(M+m)*N+n-1)...
        -A((M+m)*N+n,(M+m)*N+n+1)-A((M+m)*N+n,(M+m-1)*N+n)...
        -A((M+m)*N+n,(M+m+1)*N+n)+K2*k-K3*k*u(m*N+n)/(K4+u(m*N+n));
    end
  end
return