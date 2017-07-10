function [ P ] = EM_EStep(K,d,n,W, Mus,Sigma,score )
% Computers P matrix which is the probability that each observation is in a
% given category.

 c = (2 * pi)^(.5*d);
    S = zeros(1,K);
    inv_sig = zeros(d,d,K);
    for i = 1:K
        A = Sigma(:,:,i);
        b = eye(size(Sigma,2));
        inv_sig(:,:,i) = A\b;

         S(i) = sqrt(det(Sigma(:,:,i)));

    end
    
 

     
       P = zeros(n,K); 
for i=1:n    

    for j=1:K  
        dXM = score(i,1:d)' - Mus(j,:)';
        pl = exp(-0.5*dXM'*inv_sig(:,:,j)*dXM)/(c*S(j));
        P(i,j) = W(j)*pl;
    end
        
    P(i,:) = P(i,:)/sum(P(i,:));
end  
end



