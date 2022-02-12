function [ U ] = WAVE( X,Time,TYPE,Alpha,Beta,dx,N,C )

Lambda = Beta/dx;

if (TYPE == 1)
    for i = 1:N+1
        U(i) = sin(Lambda*(X(i)-C*Time));
    end
elseif (TYPE == 2)
    for i = 1:N+1
        if((X(i)-C*Time) <= 1.0 && (X(i)-C*Time) >= 0.0 )
            U(i) = 1.0;
        else
            U(i) = 0.0;
        end
    end
elseif (TYPE == 3)
    for i = 1:N+1
        U(i) = exp(-Alpha*(((X(i)-C*Time)-1.0)^2.0))*sin(Lambda*(X(i)-C*Time));
    end
end

end

