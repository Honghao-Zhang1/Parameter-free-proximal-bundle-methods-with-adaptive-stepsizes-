function   method_2_gpb(K,lam,x_0,eps,n,A,b)
%METHOD_2_GPB GPB.
%
%   method_2_gpb(K, lam, x_0, eps, n, A, b) runs the constant-stepsize GPB
%   method on
%
%       min ||A*x-b||_1  subject to x >= 0.
%
%   This is the GPB baseline for the nonnegative l1 feasibility experiment.
%   It uses the two-cut bundle update with fixed lam.

    x = x_0;
    y = x_0;
    % counters = zeros(K,1);
    total = 1;  
    tic
    g =  A'*sign(A*x'-b);
    fbest = norm(A*y'-b,1);
    fx = fbest;
    for k = 1:K 
        center = x; 
        a1 = g;
        b1 = fx - a1'*center'; 
        x = max(center - lam*g',0);
        % a2=zeros(n,1);
        % b2 = -Inf;
        % [~,x] =Twocut_subproblem1(n,lam,a1,a2,b1,b2,center);
        g =  A'*sign(A*x'-b);
        fx = norm(A*x'-b,1);
        if fx< eps
            cpu = toc;
            disp(['k=',num2str(k), ' total = ',num2str(total), ' cpu=', num2str(cpu),' lam=', num2str(lam), ' avginn = ', num2str(total/k)])
            % counters(k) = j;
            return
        end
        dist = norm(x-center)^2;
        Gamma =  a1'*x'+b1 ;
        m = Gamma+ 1/(2*lam)*dist;  
        if  fx < fbest
            y = x;  
            fbest = fx;
        end
        t = fbest- m;
        % t0 = t;
        % Gamma = max{a1x+b1,a2x+b2} ; Gamma1 = max{lm(;,x0),lm(;x1)}
        a2 = g ;
        b2 = fx - a2'*x';
        total = total +1;   
        j = 1;
                            time= toc;
            if time > 3600
                disp(['raito = ', num2str(fbest/eps)])
                disp('----out of time----')
                return
            end
        while t  >   eps/2
            j = j+1; 
            if (a1-a2)'*max(center -lam*a2',0)'+b1-b2 <=0
                alpha = 0;
            elseif (a1-a2)'*max(center -lam*a1',0)'+b1-b2 >=0
                alpha = 1;
            else
                left = 0;
                right = 1;
            % error = 1; 
                iter = 0;
            %bisection
                while right - left > 1e-6 && iter < 20            
                    alpha = (left+right)/2;
                    iter = iter+1;
                    error =(a1-a2)'*max(center -lam*(alpha *a1' + (1-alpha)*a2'),0)'+b1-b2;
                    if error < 0
                        right = alpha;
                    else
                        left = alpha;
                    end               
                end 
           end
      
            x = max(center -lam*(alpha *a1' + (1-alpha)*a2'),0);
            % [alpha,x] =Twocut_subproblem1(n,lam,a1,a2,b1,b2,center);
            fx = norm(A*x'-b,1);
            g =  A'*sign(A*x'-b);
            if fx < eps
                cpu = toc;
                disp(['k=',num2str(k), ' total = ',num2str(total), ' cpu=', num2str(cpu),' lam=', num2str(lam), 'avginn = ', num2str(total/k)])
                % counters(k) = j;
                return 
            end
            dist = norm(x-center)^2;
            Gamma = max(a1'*x'+b1,a2'*x'+b2);
            m = Gamma+ 1/(2*lam)*dist; 
            if  fx < fbest          
                y = x;  
                fbest = fx;
            end
            t =  fbest- m;
            % update bundle
            % Gamma = max{a1x+b1,a2x+b2}
            a1 = alpha*a1 + (1-alpha)*a2;
            a2 = g ; 
            b1 = alpha*b1 + (1-alpha)*b2;
            b2 = fx - a2'*x';
            total = total +1;
                    time = toc;
            if time > 3600
                disp(['raito = ', num2str(fbest/eps)])
                disp('----out of time----')
                return
            end
        end 

    end
    total
    k
    cpu = toc
end
