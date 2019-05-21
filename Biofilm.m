function [X,Y,S,U]=Biofilm(T,hx,hy,k,d1,d2,K1,K2,K3,K4,Neu)
  [t,u1,X,Y,M,N,K,Rx2,Ry2,A]=Constants(T,hx,hy,k,d1,d2,Neu);
  for i=1:K
    u1=new_approximation(u1,k,M,N,Rx2,Ry2,A,K1,K2,K3,K4,Neu);
  end
  U=reshape(u1(M*N+1:2*M*N),N,M); 
  S=reshape(u1(1:M*N),N,M); 
  PlotSolution(X,Y,U,t(i));  
  PlotSolution(X,Y,S,t(i));
return