function [ d ] = EM_sqdist( Mus, b )
% Computes the Euclidean distance between all points and the means

if size(Mus,1) == 1
    d = repmat(Mus',1,length(b)) - repmat(b,length(Mus),1); 
    d = d.^2;
else
    aa = sum(Mus.*Mus); bb = sum(b.*b); 
    ab = (Mus'*b);
    d = abs(repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab);
end

end

