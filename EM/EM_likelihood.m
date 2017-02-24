function [ L ] = EM_likelihood( score, W, Mus, Sigma )
% Calculates the likelihood of the current solution

[K,f] = size(Mus);
[n,d] = size(score(:,1:f));
L = 0;
U = mean(Mus);
Score_cov = cov(score(:,1:d));

for i=1:K
    A = Sigma(:,:,i);
    b = eye(size(Sigma,2)); 
    iV = A\b; 

    L = L + W(i)*(-0.5*n*log(det(2*pi*Sigma(:,:,1))) ...
        -0.5*(n-1)*(trace(iV*Score_cov)+(U-Mus(i,:))*iV*(U-Mus(i,:))'));
end



end 

    