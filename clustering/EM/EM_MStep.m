function [ W,Mus,Sigma ] = EM_MStep(n,d,K,P,score)
% Re-computes W, Mus, and Sigma values based on the updated P matrix

   Psum = sum(P,1);
   Sigma = zeros(d,d,K);
   Mus = zeros(K,d);
   W = zeros(K,1);
   
 for i = 1:K
        
%     Update weights
    
    W(i) = W(i) + Psum(i) / n;
    
    % Update Means
    
    Mus(i,:) = Mus(i,:) +  P(:,i)' * score(:,1:d) / Psum(i);
   
    % Update Covariance Matrices

%     Sigma(:,:,i) = (score(:,1:d) - repmat(Mus(i,:),n,1))' .* ...
%         repmat(P(:,i),1,d)' * (score(:,1:d) - repmat(Mus(i,:),n,1)) ...
%         ./ repmat(Psum(i),d,d);
%     

  
    for j=1:n
        dXM = score(j,1:d)'-Mus(i,:)';
        Sigma(:,:,i) = Sigma(:,:,i) + P(j,i) * dXM * dXM';
    end
    Sigma(:,:,i) = Sigma(:,:,i)/Psum(i);

%     d_XM = (score(:,1:d)' - repmat(Mus(i,:),n,1))';
%     Sigma(:,:,i) = Sigma(:,:,i) + Psum(i) * (d_XM' * d_XM);
%     Sigma(:,:,i) = Sigma(:,:,i) ./ repmat(Psum(i),d,d);

    
%  Code for initializing    
%    b = find(I == i);
%    if length(b) > d*2;
%        S = cov(score(b,1:d));
%    else 
%        S = cov(score(b,1:d));
%    end
%    Sigma(:,:,i) = S;

    
end

end

