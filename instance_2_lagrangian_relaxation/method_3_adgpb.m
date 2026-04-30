function [pi,fmax] = method_3_adgpb(q,W,T,h,x,subopt,pi_0,tau)
    % tau = 0.95;
    opt = -subopt;
    K =  300;
    bad= 0;
    pi = pi_0;
    % bad_old = 0;
    % xopt = zeros(n,K);
    y = pi;
    % counters = zeros(K,1);
    total = 0;  
    % tic
    [fmax,g] = helper_dual_oracle(q,W,T,h,pi,x);
    fbest = -fmax;
     eps = 1e-4;
     % fbest-opt
    if fbest-opt< eps
            % cpu = toc;
            % fmax = -fbest;
            % disp(['k=',num2str(0), ' total = ',num2str(total),', cpu=', num2str(cpu)])
            fmax = -fbest;
            % counters(k) = j;
            return
    end
    
    % fbest = norm(A*y'-b,1);
    fpi = fbest;
    lam = (abs(fpi-opt))/norm(g)^2;
    for k = 1:K 
                    % lam =40*(fpi-opt)/norm(g)^2;
        center = pi; 
        a1 = -g;
        b1 = fpi - a1'*center; 
        pi = center - lam*a1;
        [fmax,g] = helper_dual_oracle(q,W,T,h,pi,x);
        fpi = -fmax;
        total = total +1;
        if fpi-opt< eps
       
            % cpu = toc;
             fmax = -fbest;
            % disp(['k=',num2str(k), ' total = ',num2str(total),'cpu=', num2str(cpu)])
            % counters(k) = j;
            return
        end
        dist = norm(pi-center)^2;
        Gamma =  a1'*pi+b1 ;
        m = Gamma+ 1/(2*lam)*dist;  
        if  fpi < fbest
            y = pi;  
            fbest = fpi;
        end
        t = fbest- m;
        t_old = t;
        % t0 = t;
        % Gamma = max{a1x+b1,a2x+b2} ; Gamma1 = max{lm(;,x0),lm(;x1)}
        a2 = -g ;
        b2 = fpi - a2'*pi;
 
        j = 1;
        % eps/4+fbest/2
        while t  >   eps/4+(fbest-opt)/2          
            total = total +1;
            j = j+1; 
            % left = 0;
            % right = 1;
            % % error = 1; 
            % iter = 0;
            % % bisection
            % while right - left > 1e-6 && iter < 15            
            %     alpha = (left+right)/2;
            %     iter = iter+1;
            %     error =(a1-a2)'*(center -lam*(alpha *a1 + (1-alpha)*a2))+b1-b2;
            %     if error*((a1-a2)'* (center -lam*(left *a1 + (1-left)*a2))+b1-b2) < 0
            %         right = alpha;
            %     else
            %         left = alpha;
            %     end               
            % end     
            alpha =  ((a1 - a2)'*center- lam*a1'*a2+lam *norm(a2)^2+b1-b2)/(lam*norm(a1-a2)^2);
            if alpha > 1
                alpha = 1;
            elseif alpha < 0
                alpha = 0;
            end
            % (a1-a2)'*max(center -lam*(alpha *a1' + (1-alpha)*a2'),0)'+b1-b2
            pi= center - lam*(alpha *a1 + (1-alpha)*a2);
            [fmax,g] = helper_dual_oracle(q,W,T,h,pi,x);
            fpi = -fmax;
            if fpi-opt < eps
                % cpu = toc;
                fmax = -fbest;
                % counters(k) = j;
                % disp(['k=',num2str(k), ' total = ',num2str(total), ', cpu=', num2str(cpu)])
                return 
            end
            dist = norm(pi-center)^2;
            Gamma = max(a1'*pi+b1,a2'*pi+b2);
            m = Gamma+ 1/(2*lam)*dist; 
            if  fpi < fbest          
                y = pi;  
                fbest = fpi;
            end
            t =  fbest- m;
            if t- tau*t_old > (1-tau)*(eps/4+(fbest-opt)/4)          
                bad=1;
                lam = lam/2;
            end
            t_old = t;
            % update bundle
            % Gamma = max{a1x+b1,a2x+b2}
            a1 = alpha*a1 + (1-alpha)*a2;
            a2 = -g ; 
            b1 = alpha*b1 + (1-alpha)*b2;
            b2 = fpi - a2'*pi;
        end 
        if bad == 0  
           lam = 2*lam;
        end
        % if  bad_old ==0 && bad ==1
        %     k
        % end
        % bad_old = bad;
        % counters(k) = j;
        % ra = barm/m; 
        % xopt(:,k) = y;
    end
    % total
    % k
    % cpu = toc
end