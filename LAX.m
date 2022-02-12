function [ u,G ] = LAX( u,C,n,N,dt,dx,NTime,Sigma,Beta )

while( n <= NTime )
    for i = 2:N
        u(n+1,i) = 0.5*(u(n,i+1)+u(n,i-1))-0.5*C*(u(n,i+1)-u(n,i-1));
    end
    u(n+1,1) = 0.5*(u(n,2)+u(n,N+1))-0.5*C*(u(n,2)-u(n,N+1));
    u(n+1,N+1) = 0.5*(u(n,1)+u(n,N))-0.5*C*(u(n,1)-u(n,N));
    u(n+1,N+1) = u(n,N+1)-(Sigma/2)*(u(n,1)-u(n,N))+((Sigma^2)/2)*(u(n,1)-2*u(n,N+1)+u(n,N));
    n = n+1;
end

G = cos(Beta)-1i*Sigma*sin(Beta);

end

