function   x  =  solve_Lp( y, lambda, p )

J     =   10;
tau   =  (2*lambda.*(1-p)).^(1/(2-p)) + p*lambda.*(2*(1-p)*lambda).^((p-1)/(2-p));
x     =   zeros( size(y) );
i0    =   find( abs(y) > tau );


if ~isempty(i0)
    y0    =   y(i0);
    t     =   abs(y0);
    for  j  =  1 : J
        if length(lambda) == 1
            t    =  abs(y0) - p*lambda.*(t).^(p-1);
        else
            t = abs(y0) - p*lambda(i0).*(t).^(p-1);
        end
    end
    x(i0)   =  sign(y0).*t;
end