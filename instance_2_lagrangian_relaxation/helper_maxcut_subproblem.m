function solution = helper_maxcut_subproblem(center,l,grad,fvalue,C,n)
clear prob;

 % [~,C]=size(fvalue);
prob.buc=ones(C,1);
prob.blc=[];
prob.blx=[];
prob.bux = [];

% prob.blc(C+1)=D;
prob.c = -center'; 

prob.qosubi=1:n;
prob.qosubj=1:n;
    
prob.qoval= ones(n,1);



prob.a = zeros(C,n);
for i=1:C
    prob.a(i,:) =  grad(:,i)';   
    prob.buc(i) = l-fvalue(i) ;
    % prob.buc(i) = Inf;
end
% prob.a(C+1,:) = ones(1,n);
% prob.buc(C+1)= l;
% prob.blc(C+1)= - Inf;
% prob.blx(1) = - Inf;
% prob.bux(1) = l;
% prob.bux=[];

[~,res]=mosekopt('minimize echo(0)', prob);

% res.sol.itr.xx(1)-l
solution=res.sol.itr.xx(1:n)';





end