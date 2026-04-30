function method_6_bg(x0,beta,eps,A,b,K,n,xopt)
%METHOD_6_BG GCPBBG.
%
%   method_6_bg(x0, beta, eps, A, b, K, n, xopt) runs a diagnostic
%   classical bundle variant that uses the known generated optimizer xopt to
%   choose the prox stepsize.
%
%   This is the GCPBBG routine for the nonnegative l1 feasibility
%   experiment. It uses the known generated optimizer xopt as an oracle-like
%   benchmark for choosing the prox stepsize.

    total = 1;
    % counters = zeros(K,1);
    x= x0;
    xopt = xopt(:);
    fbest = norm(A*x'-b,1);
    tic
    for k = 1:K  
        g =  A'*sign(A*x'-b);
        center = x; 
        fx = norm(A*x'-b,1);
        fstar = norm(A * xopt - b, 1);
        gap_center = max(fx - fstar, 1e-16);
        a1 = g;
        b1 = fx - a1'*center'; 
        % a2=zeros(n,1);
        % b2=-Inf;
        % Impractical rule: use the true optimizer x_* to set the prox scale
        % from the current serious iterate.
        lam = norm(center' - xopt)^2 / gap_center;
        % [~,x] = Twocut_subproblem1(n,lam,a1,a2,b1,b2,center);
         x = max(center - lam*g',0);
        g =  A'*sign(A*x'-b);
        fx = norm(A*x'-b,1);
        if  fx < fbest
            fbest = fx;
        end
        if fx< eps  
            cpu = toc;
            % counters(k) = j;
            disp(['k=',num2str(k), ' total = ',num2str(total), ' cpu=', num2str(cpu), 'avginner=', num2str(total/k)])
            return
        end
        Gamma =  a1'*x'+b1 ;
        % Gamma = max{a1x+b1,a2x+b2} ; Gamma1 = max{lm(;,x0),lm(;x1)}
        a2 = g ;
        b2 = fx - a2'*x';
        total = total +1;
        j = 1;
        % time= toc;
        % if time > 3600
        %     disp('----out of time----')
        %     return
        % end
        while beta*(gap_center - (Gamma - fstar)) > gap_center - (fx - fstar)
            total = total +1;
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
            % [alpha,x] = Twocut_subproblem1(n,lam,a1,a2,b1,b2,center);
            x = max(center -lam*(alpha *a1' + (1-alpha)*a2'),0);
            fx = norm(A*x'-b,1);
            if  fx < fbest
                    fbest = fx;
            end
            g =  A'*sign(A*x'-b);
            if fx < eps
                cpu = toc;
                % counters(k) = j;
                disp(['k=',num2str(k), ' total = ',num2str(total), ' cpu=', num2str(cpu), 'avginner=', num2str(total/k)])
                return 
            end
            Gamma = max(a1'*x'+b1,a2'*x'+b2); 
            % update bundle
            % Gamma = max{a1x+b1,a2x+b2}
            a1 = alpha*a1 + (1-alpha)*a2;
            a2 = g ; 
            b1 = alpha*b1 + (1-alpha)*b2;
            b2 = fx - a2'*x';
        end 
                    time= toc;
            if time > 3600
                disp(['raito = ', num2str(fbest/eps)])
                disp('----out of time----')
                return
            end

    end 


end
