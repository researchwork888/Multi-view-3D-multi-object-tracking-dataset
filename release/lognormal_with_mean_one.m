function [mean,std_dev]=lognormal_with_mean_one(percen)

% input : lognormal std deviation 
    mu = @(sigma) -sigma^2/2; 
    percen_v = percen^2;
    sigmaa = sqrt(log(percen_v+1));
    std_dev = sigmaa;
    mean = mu(std_dev);
    
    
    


end