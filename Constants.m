function [t,u1,X,Y,M,N,K,Rx2,Ry2,A]=Constants(T,hx,hy,k,d1,d2,Neu)
  x=0:hx:1; y=0:hy:1; t=0:k:T;
  M=length(x); N=length(y); K=length(t);
  [X,Y]=meshgrid(x,y);
  Rx1=d1*k/hx^2; Ry1=d1*k/hy^2;
  Rx2=d2*k/hx^2; Ry2=d2*k/hy^2;
  A=constant_matrix(M,N,Rx1,Ry1,Neu);
  u1=initial_profile(X,Y,M,N);
return