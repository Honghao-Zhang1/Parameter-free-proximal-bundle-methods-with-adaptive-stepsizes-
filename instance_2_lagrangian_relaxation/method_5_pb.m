  function method_5_pb(pi_0,beta,eps,lammin,a,q,W,T,h,x,opt,K)
    total = 1;
    opt = -opt;
    % counters = zeros(K,1);
    pi= pi_0;
    index = 0;
    % tic
    [fmax,g] = helper_dual_oracle(q,W,T,h,pi,x);
    fbest = -fmax;
    fvalue = fbest;
    lam = (abs(fvalue-opt))/norm(g)^2;
    for k = 1:K  
        center = pi; 
        a1 = -g;
        fcenter = fvalue;
        b1 = fcenter - a1'*center; 
        % a2=zeros(n,1);
        % b2=-Inf;
        % [~,x] = Twocut_subproblem1(n,lam,a1,a2,b1,b2,center);
         pi = center - lam*a1;
         [fmax,g] = helper_dual_oracle(q,W,T,h,pi,x);
         fvalue = -fmax;
        if  fvalue < fbest
           fbest = fvalue;
        end
        if fbest-opt< eps  
            % cpu = toc;
            % counters(k) = j;
            % disp(['k=',num2str(k), ' total = ',num2str(total), ' cpu=', num2str(cpu), 'avginner=', num2str(total/k)])
            return
        end
        Gamma =  a1'*pi+b1 ;
        % Gamma = max{a1x+b1,a2x+b2} ; Gamma1 = max{lm(;,x0),lm(;x1)}
        a2 = -g ;
        b2 = fvalue - a2'*pi;
        total = total +1;
        j = 1;
        % time= toc;
        % if time > 1800
        %       disp('----out of time----')
        %     return
        % end
        while beta*(fcenter -Gamma )  >   fcenter - fvalue
            lam_aux = 2*lam *(1+(fcenter-fvalue)/(fcenter-Gamma));
            lam = min(lam,max([lam_aux,lam/a,lammin]));
            total = total +1;
            j = j+1; 

            alpha =  ((a1 - a2)'*center- lam*a1'*a2+lam *norm(a2)^2+b1-b2)/(lam*norm(a1-a2)^2);
            if alpha > 1
                alpha = 1;
            elseif alpha < 0
                alpha = 0;
            end
    
            pi= center - lam*(alpha *a1 + (1-alpha)*a2);
            [fmax,g] = helper_dual_oracle(q,W,T,h,pi,x);
            fvalue = -fmax;
            if  fvalue < fbest
                    fbest = fvalue;
            end
            if fbest-opt < eps
                % cpu = toc;
                % counters(k) = j;
                % disp(['k=',num2str(k), ' total = ',num2str(total), ' cpu=', num2str(cpu), 'avginner=', num2str(sum(total)/k)])
                return 
            end
            Gamma = max(a1'*pi+b1,a2'*pi+b2); 
            % update bundle
            % Gamma = max{a1x+b1,a2x+b2}
            a1 = alpha*a1 + (1-alpha)*a2;
            a2 = -g ; 
            b1 = alpha*b1 + (1-alpha)*b2;
            b2 = fvalue - a2'*pi;
            % time= toc;
            % if time > 1800
            %     disp(['raito = ', num2str((fbest-opt)/eps)])
            %     disp('----out of time----')
            %     return
            % end
        end 
        if j==1
            index= index+1;
        end
        if index >= 5
            lam_aux = a*lam_aux;
        else
            lam_aux = 2*lam *(1+(fcenter-fvalue)/(fcenter-Gamma));
        end
        lam = min(lam_aux,10*lam);
        % counters(k) = j;
        % ra = barm/m; 
        % xopt(:,k) = y;
    end 


end