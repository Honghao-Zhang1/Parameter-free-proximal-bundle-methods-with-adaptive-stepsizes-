function  method_8_adgpb(K,lam,x_0,eps,tau,n,A,b,D)
    bad= 0;
    % bad_old = 0;
    % x = zeros(n,K);
    x = x_0;
    y = x_0;
    % counters = zeros(K,1);
    total = 1;  
    tic
    g =  A'*sign(A*x'-b);
    fbest = norm(A*y'-b,1);
    fx = fbest;
    if all(g>=0)
        Lb = fx-g'*x' ;
    else
        Lb = fx-g'*x' + D*min(g);
    end
    finitial = fbest;
    Lbini = Lb;
    Affine_g = zeros(1,n);
    Affine_fvalue = zeros(1);
    phik=zeros(1);
    hatgk = zeros(1);
    lamsum = 0;
    beta = 1/2;
    lamini = lam;
    lammin = lam;
    lammax = lam;
    for k = 1:K 
        center = x; 
        a1 = g;
        b1 = fx - a1'*center'; 
        fcenter = fbest;
        % a2=zeros(n,1);
        % b2 = -Inf;
        % [~,x] = Twocut_subproblem1(n,lam,a1,a2,b1,b2,center);
        x = helper_project_simplex(center'- lam*a1,D)';
         g =  A'*sign(A*x'-b);
        fx = norm(A*x'-b,1);
        if fx-Lb< eps
            cpu = toc;
            disp(['k=',num2str(k), ', total = ',num2str(total), ', cpu=', num2str(cpu), ', avginner= ', num2str(total/k),  ', betaratio = ', num2str(2*beta)])
            disp([', lamratio =', num2str(lam/lamini),', minlamratio =', num2str(lammin/lammin), ', maxlamratio =', num2str(lammax/lamini)])
            %disp(['*PG =', num2str(fbest/(finitial-Lbini)), ', DG=', num2str((-Lb)/(finitial-Lbini)),' PDG =', num2str((fbest-Lb)/(finitial-Lbini)), 'DG/PG=', num2str((-Lb)/fbest)])
            % counters(k) = j;
            disp(['fx=',num2str(fx), ', e^Tx=', num2str(ones(1,n)*x')])
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
        while t  >   eps/4+beta*(fcenter-Lb)
            j = j+1; 
            % left = 0;
            % right = 1;
            % % error = 1; 
            % iter = 0;
            %bisection
            if (a1-a2)'*helper_project_simplex(center' -lam*a2,D)+b1-b2 <=0
                alpha = 0;
            elseif (a1-a2)'*helper_project_simplex(center' -lam*a1,D)+b1-b2 >=0
                alpha = 1;
            else
                left = 0;
                right = 1;
            % error = 1; 
                iter = 0;
                while right - left > 1e-6 && iter < 20            
                    alpha = (left+right)/2;
                    iter = iter+1;
                    error =(a1-a2)'*helper_project_simplex(center' -lam*(alpha *a1 + (1-alpha)*a2),D)+b1-b2;
                    if error < 0
                        right = alpha;
                    else
                        left = alpha;
                    end               
                end  
            end
            % alpha =  ((a1 - a2)'*center'- lam*a1'*a2+lam *norm(a2)^2+b1-b2)/(lam*norm(a1-a2)^2);
            % if alpha > 1
            %     alpha = 1;
            % elseif alpha < 0
            %     alpha = 0;
            % end
            % (a1-a2)'*max(center -lam*(alpha *a1' + (1-alpha)*a2'),0)'+b1-b2
            % [alpha,x] =Twocut_subproblem1(n,lam,a1,a2,b1,b2,center);
            x = helper_project_simplex(center' -lam*(alpha *a1 + (1-alpha)*a2),D)';
            fx = norm(A*x'-b,1);
            g =  A'*sign(A*x'-b);
            if fx-Lb < eps
                cpu = toc;
                fx
                % counters(k) = j;
                disp(['k=',num2str(k), ', total = ',num2str(total), ', cpu=', num2str(cpu), ', avginner= ', num2str(total/k), ', lamratio =', num2str(lam/lamini), ', betaratio = ', num2str(2*beta),', minlamratio =', num2str(lammin/lammin), ', maxlamratio =', num2str(lammax/lamini)])
                %disp(['*PG =', num2str(fbest/(finitial-Lbini)), ', DG=', num2str((-Lb)/(finitial-Lbini)),' PDG =', num2str((fbest-Lb)/(finitial-Lbini)), 'DG/PG=', num2str((-Lb)/fbest)])
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
            if t- tau*t_old > (1-tau)*(eps/8+beta*(fcenter-Lb)/2)           
                 bad=1;
                lam = lam/2;
                if lam < lammin
                    lammin = lam;
                end
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
        if time > 1800
            %disp(['k=',num2str(k), ', total = ',num2str(total), ', cpu=', num2str(time), ', avginner= ', num2str(total/k), ', lamratio =', num2str(lam/lamini), ', betaratio = ', num2str(2*beta),', minlamratio =', num2str(lammin/lammin), ', maxlamratio =', num2str(lammax/lamini)])
            disp(['k=',num2str(k), ', total = ',num2str(total), ', cpu=', num2str(time), ', avginner= ', num2str(total/k),  ', betaratio = ', num2str(2*beta)])
            disp([', lamratio =', num2str(lam/lamini),', minlamratio =', num2str(lammin/lammin), ', maxlamratio =', num2str(lammax/lamini)])
            %disp(['*PG =', num2str(fbest/(finitial-Lbini)), ', DG=', num2str((-Lb)/(finitial-Lbini)),' PDG =', num2str((fbest-Lb)/(finitial-Lbini)), 'DG/PG=', num2str((-Lb)/fbest)])
            % counters(k) = j;
            disp(['fx=',num2str(fbest), ', e^Tx=', num2str(ones(1,n)*x')])
            disp(['*PG =', num2str(fbest/(finitial-Lbini)), ', DG=', num2str((-Lb)/(finitial-Lbini)),' PDG =', num2str((fbest-Lb)/(finitial-Lbini)), 'DG/PG=', num2str((-Lb)/fbest)])
            return
        end
        LBg = 1/lam*(center-x);
        Lbvalue = Gamma-1/lam*(center-x)*x';
        lamsum = lamsum+lam;
 
        Affine_g = Affine_g+lam*LBg;
        Affine_fvalue = Affine_fvalue +lam*Lbvalue;
        phik = phik + lam*fbest;
        hatgk = hatgk + beta*lam*(fbest-Lb);
        if checkPowerOfTwo(k) 
            hatgk(end+1) =  beta*lam*(fbest-Lb);
            lamsum(end+1) = lam;
            phik(end+1) = lam*fbest;
            Affine_g(end+1,:) = lam* LBg;
            Affine_fvalue(end+1)= lam*Lbvalue;
        end
        for j=1:length(lamsum)
            % lamsum(j)
            % Affine_fvalue(j)/lamsum(j)
            % con = Affine_fvalue
            % normg = norm(Affine_g)
            if all(Affine_g(j,:)>=0)
                    NLB  = Affine_fvalue(j)/lamsum(j);
            else
                   NLB =Affine_fvalue(j)/lamsum(j)+D*min(Affine_g(j,:)/lamsum(j));
            end
            % NLB =Affine_fvalue(j)/lamsum(j)+D*min(Affine_g(j,:)/lamsum(j));
            Lb = max(NLB,Lb);
            if (j == length(lamsum)-1) &  (2* hatgk> phik(j) - Lb*lamsum(j))
                beta = beta/2;
            end
        end
        if bad == 0  
           lam = 2*lam;
        end
        if lam >lammax
            lammax = lam;
        end
        % if  bad_old ==0 && bad ==1
        %     k
        % end
        % bad_old = bad;
        % counters(k) = j;
        % ra = barm/m; 
        % x(:,k) = y;
    end
    % total
    % k
    % cpu = toc
end

   function isPowerOfTwo = checkPowerOfTwo(n)
    if n >1 && mod(n, 2) == 0
        isPowerOfTwo = (bitand(n, n - 1) == 0);
    else
        isPowerOfTwo = false;
    end
   end
