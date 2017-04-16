function [ Theta, obj ] = solve_step_1_nesterov(cx, cy, SX, SZX, LX, LZX, Lambda, mu, tol, theta)

Sx	= cx'*cx;
Sy 	= cy'*cy;
Sxy	= cx'*cy;
N 	= size(cx, 1);
q = size(cx, 2);
p = size(cy, 2);

maxiter = 1000;
eta     = 1.5;
nobj	= 10;
bconv	= 0;
obj	= zeros(maxiter, 1);
%theta	= struct;
%theta.yy = eye(p);
%theta.xy = zeros(q,p);
L	= 1;
thk_0 	= 2/3;
ls_maxiter = 300;

[obj1, init_flag]  = compute_gradient( theta, Sx, Sxy, Sy,N, 'n', SX, SZX, LX, LZX, Lambda, mu);
if init_flag == 1 && verbose == true
    fprintf('sCGGM: error! initial Theta_yy not positive definite!\n');
end
obj(1) = obj1;
xk      = theta;
zk   	= theta;
thk     = thk_0;

for iter = 2:maxiter
    thk  = (sqrt( thk^4 + 4 * thk^2 ) - thk^2) / 2; % momentum of acceleration
    y.xy = (1 - thk) * xk.xy + thk * zk.xy;
    y.yy = (1 - thk) * xk.yy + thk * zk.yy;
    [ fyk, flagy, grady] = compute_gradient( y, Sx, Sxy, Sy, N, 'y', SX, SZX, LX, LZX, Lambda, mu); % compute the objective and gradient for y
    
    % line search
    ik = 0;
    while ( 1 )
        % gradient step
        zk_grady.xy = zk.xy - 1/(L*thk) * grady.xy;
        zk_grady.yy = zk.yy - 1/(L*thk) * grady.yy;
        % proximal step
        %zk1 		= compute_prox_l1( zk_grady, 2*lambda1/(thk*L)) ;
        zk1 = zk_grady;
        % gradient step
        y_grady.xy	= y.xy - 1/L * grady.xy;
        y_grady.yy	= y.yy - 1/L * grady.yy;
        % proximal step
        %xk1         = compute_prox_l1( y_grady, 2*lambda1/(L));
        xk1 = y_grady;
        [fxk1, flagxk1] = compute_gradient(xk1, Sx, Sxy, Sy, N ,'n', SX, SZX, LX, LZX, Lambda, mu);
        [~, flagzk1]    = chol(zk1.yy);
        
        if ( flagzk1 == 0 && flagy ==0 && flagxk1 ==0 ) % xk1,zk1,y all positive definite
            xk1_y.xy    = xk1.xy - y.xy;
            xk1_y.yy    = xk1.yy - y.yy;
            lfxk1_y     = fyk + grady.xy(:)'* (xk1_y.xy(:)) + grady.yy(:)'*(xk1_y.yy(:));
            diffxk1y.xy = xk1.xy - y.xy;
            diffxk1y.yy = xk1.yy - y.yy;
            RHS         = lfxk1_y + L/2 *(sum(diffxk1y.xy(:).^2) + sum(diffxk1y.yy(:).^2));
            if fxk1 <= RHS + tol
                xk = xk1;
                zk = zk1;
                bconv = 1;
                break; % line search converged
            end
        end
        
        ik = ik + 1;
        
        if ( ik > ls_maxiter )
            if verbose
                fprintf( 'sCGGM: line search not converging,ik = %d\n',ik);
            end
            bconv = 0;
            iter  = max(1, iter - 1);
            Theta = xk;
            break;
        end
        L = L * eta;
    end
    obj(iter)  = fxk1;
    if bconv == 0
        break;
    end
    
    if ( iter > nobj + 1)
        value           = obj(iter);
        prevVals        = obj(iter - nobj);
        avgimprovement  = abs( prevVals - value )/nobj;
        relAvgImpr      = avgimprovement / abs( value ) ; % relative average improvement
        
        if ( relAvgImpr < tol )
            bconv = 1;
            break;
        end
    end
end

Theta = xk;
obj   = obj(1:iter);

end

function [ value, flag, grad] = compute_gradient(theta, Sx, Sxy, Sy, N, gradient, SX, SZX, LX, LZX, Lambda, mu)

flag        = 0;
[ cyy, p ]  = chol(theta.yy);

if ( p > 0 )
    if strcmp(gradient, 'y') == 1 && verbose
        fprintf( 'sCGGM: Theta_yy not positive definite!\n' );
    end
    flag        = 1;
    value       = inf;
    grad        = theta;
    return;
end

logdetyy = 2 * sum(log(diag(cyy) ));

if ( isnan(logdetyy) || isinf(logdetyy) )
    if verbose
        fprintf( 'sCGGM: logdet Theta_yy is Nan or Inf!\n' );
    end
    flag = 1;
    value = inf;
    grad = theta;
    return;
end

icyy	 = cyy \ eye(size(cyy,2));
ithetayy = icyy * icyy';
txyityy  = (theta.xy)*ithetayy;
XtXth    = Sx*txyityy;
txyXtXth = (theta.xy)'*Sx*txyityy;

l1 = trace( (theta.yy)*Sy );
l2 = trace( Sxy*(theta.xy)' );
l3 = trace( txyXtXth );
value = 0.5*l1 + l2 + 0.5*l3 - 0.5*N*logdetyy ;
value = value / N;
value = value + 0.5 * mu * norm([theta.yy; theta.xy] - [SX; SZX] + [LX; LZX] + Lambda/mu, 'fro')^2;

LambdaX = Lambda(1:size(SX,1),:);
LambdaZX = Lambda((size(SX,1)+1):end,:);
if strcmp('y',gradient) ==1
    grad.xy = (Sxy + XtXth)/N + LambdaZX + mu * theta.xy - mu * SZX + mu * LZX;
    grad.yy = 0.5*(Sy - N*ithetayy - ithetayy*txyXtXth)/N + LambdaX + mu * theta.yy - mu * SX + mu * LX;
else
    grad = [];
end
end
