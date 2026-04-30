  function method_5_pb(lam,x0,beta,eps,A,b,K,n,lammin,a)
%METHOD_5_PB*.
%
%   method_5_pb(lam, x0, beta, eps, A, b, K, n, lammin, a) runs a classical
%   bundle-style method on the l1 feasibility model
%
%       min ||A*x-b||_1  subject to x >= 0.
%
%   This is the PB* routine for the nonnegative l1 feasibility experiment.
%   It uses the two-cut model with a classical sufficient decrease
%   serious/null decision.

    total = 1;
    % counters = zeros(K,1);
    x= x0;
    index = 0;
    fbest = norm(A*x'-b,1);
    tic
    for k = 1:K  
        g =  A'*sign(A*x'-b);
        center = x; 
        a1 = g;
        fcenter = norm(A*x'-b,1);
        b1 = fcenter - a1'*center'; 
        % a2=zeros(n,1);
        % b2=-Inf;
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
                          time= toc;
            if time > 3600
                disp(['raito = ', num2str(fbest/eps)])
                disp('----out of time----')
                return
            end
        while beta*(fcenter -Gamma )  >   fcenter - norm(A*x'-b,1)
            lam_aux = 2*lam *(1+(fcenter-fx)/(fcenter-Gamma));
            lam = min(lam,max([lam_aux,lam/a,lammin]));
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
                disp(['k=',num2str(k), ' total = ',num2str(total), ' cpu=', num2str(cpu), 'avginner=', num2str(sum(total)/k)])
                return 
            end
            Gamma = max(a1'*x'+b1,a2'*x'+b2); 
            % update bundle
            % Gamma = max{a1x+b1,a2x+b2}
            a1 = alpha*a1 + (1-alpha)*a2;
            a2 = g ; 
            b1 = alpha*b1 + (1-alpha)*b2;
            b2 = fx - a2'*x';
                              time= toc;
            if time > 3600
                disp(['raito = ', num2str(fbest/eps)])
                disp('----out of time----')
                return
            end
        end 

        if j==1
            index= index+1;
        end
        if index >= 5
            lam_aux = a*lam_aux;
        else
            lam_aux = 2*lam *(1+(fcenter-fx)/(fcenter-Gamma));
        end
        lam = min(lam_aux,10*lam);
        % counters(k) = j;
        % ra = barm/m; 
        % xopt(:,k) = y;
    end 


end
