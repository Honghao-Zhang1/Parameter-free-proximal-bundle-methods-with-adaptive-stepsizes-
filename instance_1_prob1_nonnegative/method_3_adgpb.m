function  method_3_adgpb(K,lam,x_0,eps,tau,n,A,b)
%METHOD_3_ADGPB Ad-GPB*.
%
%   method_3_adgpb(K, lam, x_0, eps, tau, n, A, b) applies the adaptive
%   proximal bundle method to
%
%       min ||A*x-b||_1  subject to x >= 0.
%
%   This is the Ad-GPB* routine for the nonnegative l1 feasibility
%   experiment. It uses the two-cut bundle update, adaptive null-step
%   contraction test, and end-of-cycle increase of lam after good cycles.

    bad= 0;
    initial=lam;
    % bad_old = 0;
    % xopt = zeros(n,K);
    x = x_0;
    y = x_0;
    counters = zeros(K,1);
    total = 0;  
    tic
    g =  A'*sign(A*x'-b);
    fbest = norm(A*y'-b,1);
    fx = fbest;
    for k = 1:K 
          fcenter = fbest;      
        center = x; 
        a1 = g;
        b1 = fx - a1'*center'; 
        % a2=zeros(n,1);
        % b2 = -Inf;
        % [~,x] = Twocut_subproblem1(n,lam,a1,a2,b1,b2,center);
        x = max(center- lam*a1',0);
         g =  A'*sign(A*x'-b);
        fx = norm(A*x'-b,1);
        if fx< eps
            cpu = toc;
            disp(['k=',num2str(k), ' total = ',num2str(total), ' cpu=', num2str(cpu), ' avginner= ', num2str(total/k)])
            counters(k) = j;
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
        t_old = t;
        % t0 = t;
        % Gamma = max{a1x+b1,a2x+b2} ; Gamma1 = max{lm(;,x0),lm(;x1)}
        a2 = g ;
        b2 = fx - a2'*x';
        total = total +1;
        j = 1;
        % eps/4+fbest/2
        while t  >   eps/4+fcenter/2
            
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
            fx = norm(A*x'-b,1);
            g =  A'*sign(A*x'-b);
            if fx < eps
                cpu = toc;
                % counters(k) = j;
                disp(['k=',num2str(k), ' total = ',num2str(total), ' cpu=', num2str(cpu), ' avginner= ', num2str(total/k)])
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
            % eps/8+fbest/4
            if t- tau*t_old > (1-tau)*(eps/8+fcenter/4)           
                bad=1;
                lam = lam/2;
            end
            t_old = t;
            % update bundle
            % Gamma = max{a1x+b1,a2x+b2}
            a1 = alpha*a1 + (1-alpha)*a2;
            a2 = g ; 
            b1 = alpha*b1 + (1-alpha)*b2;
            b2 = fx - a2'*x';
            total = total +1;
        end 
        time = toc;
        if time > 3600
            disp(['raito = ', num2str(fbest/eps)])
            disp("run out of time")
            return
        end
        if bad == 0  
           lam = 2*lam;
        end
  
    end
end
