functions{
	//covariance function for main portion of the model
	matrix main_GP(
		int Nx,
		vector x,
		int Ny,
		vector y, 
		real alpha1,
		real alpha2,
		real rho1,
		real rho2,
		real rho3, 
		real FR){
					matrix[Nx, Ny] K1;
					matrix[Nx, Ny] K2;
					matrix[Nx, Ny] Sigma;
	
					//specifying random Gaussian process that governs covariance matrix
					for(i in 1:Nx){
						for (j in 1:Ny){
							K1[i,j] = alpha1*exp(-square(x[i]-y[j])/2/square(rho1));
						}
					}
					
					//specifying random Gaussian process incorporates heart rate
					for(i in 1:Nx){
						for(j in 1:Ny){
							K2[i, j] = alpha2*exp(-2*square(sin(pi()*fabs(x[i]-y[j])*FR))/square(rho2))*
							exp(-square(x[i]-y[j])/2/square(rho3));
						}
					}
						
					Sigma = K1+K2;
					return Sigma;
				}
	//function for posterior calculations
	vector post_pred_rng(
		real a1,
		real a2,
		real r1, 
		real r2,
		real r3,
		real FR, 
		real sn,
		int No,
		vector xo,
		int Np, 
		vector xp,
		vector yobs){
				matrix[No,No] Ko;
				matrix[Np,Np] Kp;
				matrix[No,Np] Kop;
				matrix[Np,No] Ko_inv_t;
				vector[Np] mu_p;
				matrix[Np,Np] Tau;
				matrix[Np,Np] L2;
				vector[Np] yp;
	
	//--------------------------------------------------------------------
	//Kernel Multiple GPs for observed data
	Ko = main_GP(No, xo, No, xo, a1, a2, r1, r2, r3, FR);
	for(n in 1:No) Ko[n,n] += sn;
		
	//--------------------------------------------------------------------
	//kernel for predicted data
	Kp = main_GP(Np, xp, Np, xp,  a1, a2, r1, r2,  r3, FR);
	for(n in 1:Np) Kp[n,n] += sn;
		
	//--------------------------------------------------------------------
	//kernel for observed and predicted cross 
	Kop = main_GP(No, xo, Np, xp,  a1, a2, r1, r2, r3, FR);
	
	//--------------------------------------------------------------------
	//Algorithm 2.1 of Rassmussen and Williams... 
	Ko_inv_t = Kop'/Ko;
	mu_p = Ko_inv_t*yobs;
	Tau=Kp-Ko_inv_t*Kop;
	L2 = cholesky_decompose(Tau);
	yp = mu_p + L2*rep_vector(normal_rng(0,1), Np);
	return yp;
	}
}

data { 
	int<lower=1> N1;
	int<lower=1> N2;
	vector[N1] X; 
	vector[N1] Y;
	vector[N2] Xp;
}

transformed data { 
	vector[N1] mu;
	for(n in 1:N1) mu[n] = 0;
}

parameters {
	real<lower=0> a1;
	real<lower=0> a2;
	real<lower=0> r1;
	real<lower=0> r2;
	real<lower=0> r3;
	real<lower=0> FR; 
	real<lower=0> sigma_sq;
}

model{ 
	matrix[N1,N1] Sigma;
	matrix[N1,N1] L_S;
	
	//using GP function from above 
	Sigma = main_GP(N1, X, N1, X, a1, a2, r1, r2, r3, FR);
	for(n in 1:N1) Sigma[n,n] += sigma_sq;
	
	L_S = cholesky_decompose(Sigma);
	Y ~ multi_normal_cholesky(mu, L_S);
	
	//priors for parameters
	a1 ~ normal(0,0.5);
	a2 ~ normal(0,5);
	//incorporate minimum and maximum distances - use invgamma
	r1 ~ normal(4,1);
	r2 ~ normal(0,5);
	r3 ~ normal(10,3);
	FR ~ normal(1, 0.05);	
	sigma_sq ~ normal(0,2);
}

generated quantities {
	vector[N2] Ypred = post_pred_rng(a1, a2, r1, r2, r3, FR, sigma_sq, N1, X, N2, Xp, Y);
}

