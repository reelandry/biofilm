function A=constant_matrix(M,N,Rx1,Ry1,Neu)
  A=eye(2*M*N);
  for m=1:M-2
    A(m*N+1,m*N+2)=-Neu;         A((m+1)*N,(m+1)*N-1)=-Neu;
    A((M+m)*N+1,(M+m)*N+2)=-Neu; A((M+m+1)*N,(M+m+1)*N-1)=-Neu;
    for n=2:N-1
      A(m*N+n,m*N+n)=1+2*(Rx1+Ry1);
      A(m*N+n,m*N+n-1)=-Ry1;   A(m*N+n,m*N+n+1)=-Ry1;
      A(m*N+n,(m-1)*N+n)=-Rx1; A(m*N+n,(m+1)*N+n)=-Rx1;
    end
  end
  for n=1:N
    A(n,N+n)=-Neu;           A(M*N+1-n,(M-1)*N+1-n)=-Neu;
    A(M*N+n,(M+1)*N+n)=-Neu; A(2*M*N+1-n,(2*M-1)*N+1-n)=-Neu;
  end
return