// Utilization simulation
data {
    real<lower=0> alpha; 
    real rho;
    real sigma;
    int<lower=1> N_pre; 
    int<lower=1> N_post;
    real x_pre[N_pre]; 
    real x_post[N_post];
    real<lower=0> mu_util;
    real delta_util; 
    real<lower=0> icl;
    real<lower=0> clip_amt;
}

transformed data {
    // Setting pre clip properties
    matrix[N_pre, N_pre] K_pre = cov_exp_quad(x_pre, alpha, rho); 
    vector[N_pre] mu_pre = rep_vector(0, N_pre); 
    
    // Post clip data is less autocorrelated
    matrix[N_post, N_post] K_post = cov_exp_quad(x_post, alpha - 5, rho - 2);
    vector[N_post] mu_post = rep_vector(0, N_post);
    
    // Adding sigma to diagonal 
    for (n in 1:N_pre){
        K_pre[n, n] = K_pre[n, n] + sigma;
    }
    
    // Assuming same sigma pre/post-clip
    for (n in 1:N_post){
        K_post[n, n] = K_post[n, n] + sigma;
    }
}

parameters {
    vector[N_pre] y_pre;
    vector[N_post] y_post;
}

model {
    y_pre ~ multi_normal(mu_pre, K_pre);
    y_post ~ multi_normal(mu_post, K_post);
}

generated quantities {
    vector<lower=0, upper=icl>[N_pre] util_amt_pre;
    vector<lower=0, upper=icl+clip_amt>[N_post] util_amt_post;
    real util_rate_pre[N_pre]; 
    real util_rate_post[N_post]; 

    util_amt_pre = multi_normal_rng(mu_pre, K_pre) + rep_vector(mu_util, N_pre);
    util_amt_post = multi_normal_rng(mu_post, K_post) + rep_vector(mu_util + delta_util, N_post);
    
    for (n in 1:N_pre) {
        util_rate_pre[n] = util_amt_pre[n] / icl;
    }
    
    for (n in 1:N_post) {
        util_rate_post[n] = util_amt_post[n] / (icl + clip_amt);
    }
}
