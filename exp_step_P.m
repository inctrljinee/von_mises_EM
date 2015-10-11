function [E]=exp_step_P(xt,X,k,alpha,MX,m)
%%%
%%% EM Algorithm with Von-mises Distribution by Spturtle (Taejin Park)
%%% 
xs=size(X);
E=zeros(xs(1),k);

   
for j=1:k
        
%%% Delta 
dXM = X(:)-MX(j);

%%% Expectation
E(:,j) = alpha(j)*(   1/( sum ( exp(m(j)*cos(xt) ) ) )*exp(m(j)*cos(dXM)  )    );

end

for i=1:xs(1) 
% Normalize Expectation value(Bayes Rule Denominator)
E(i,:) = E(i,:)/sum(E(i,:));
end


end