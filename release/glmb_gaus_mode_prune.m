function [log_w_new,x_new,P_new]= glmb_gaus_mode_prune(log_w,x,P,elim_threshold)
w = exp(log_w); %un-log it first !!!!

idx= find( w > elim_threshold );
w_new= w(:,idx);
x_new= x(:,idx);
P_new= P(:,:,idx);
log_w_new = log(w_new);
