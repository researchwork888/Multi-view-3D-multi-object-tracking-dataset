function [log_w_new,x_new,P_new]= glmb_gaus_mode_merge(log_w,x,P,threshold)

 
w = exp(log_w) ; % %%%%% un-log it first !!!
L= length(w); x_dim= size(x,1);
I= 1:L;
el= 1;

if all(w==0)
    w_new = [];
    x_new = [];
    P_new = [];
    return;
end

while ~isempty(I),
    [notused,j]= max(w); j= j(1);
    Ij= []; iPt= inv(P(:,:,j));
    w_new(el,1)= 0; 
    x_new(:,el)= zeros(x_dim,1); P_new(:,:,el)= zeros(x_dim,x_dim);
    for i= I
        val= (x(:,i)-x(:,j))'*iPt*(x(:,i)-x(:,j));
        if val <= threshold,
            Ij= [ Ij i ];
        end;
    end;
   
   
   w_new(el,1)= sum(w(Ij));
   x_new(:,el)= wsumvec(w(Ij),x(:,Ij),x_dim);
   P_new(:,:,el)= wsummat(w(Ij),P(:,:,Ij),x_dim);

   x_new(:,el)= x_new(:,el)/w_new(el);
   P_new(:,:,el)= P_new(:,:,el)/w_new(el);
   I= setdiff(I,Ij);
   w(Ij)= -1;
   el= el+1;
end

if size(w_new,2) == 1
    w_new = w_new';
end
log_w_new = log(w_new);


function out = wsumvec(w,vecstack,xdim)
    wmat = repmat(w,[xdim,1]);
    out  = sum(wmat.*vecstack,2);

function out = wsummat(w,matstack,xdim)
    w = reshape(w',[1,1,size(w')]);
    wmat = repmat(w,[xdim,xdim,1]);
    out = sum(wmat.*matstack,3);
    
