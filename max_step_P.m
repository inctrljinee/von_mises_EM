function  [alpha,M_new,m]=max_step_P(xt,X,k,E,M,eps,A)
%%%
%%% EM Algorithm with Von-mises Distribution by Spturtle (Taejin Park)
%%% 
alpha = zeros(1,k); 
M_new=zeros(1,k);

 m=zeros(1,k);
 
 mix=0.1;
for i=1:k,  % Compute weights

M_new(i)=mix * M(i) + (1-mix)*atan( sum(  sin (X(:,1)).*E(:,i)  )/ sum(  cos(X(:,1)).*E(:,i)  ) );
val=sum(  cos(X(:,1)-M_new(i) ).*E(:,i)  )/sum( E(:,i) ); 
B=(A-val);
[c,ind]=min( abs(B) );
m(i)=-100 + 0.01*(min(ind)-1);
alpha(i) = mean(E(:,i));

end

end