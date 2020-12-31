functions{
	//covariance function for main portion of the model
	matrix main_GP(
		int Nx,
		vector x,
		int Ny,
		vector y, 
		real alpha1,
		real alpha2,
		real alpha3,
		real rho1,
		real rho2,
		real rho3, 
		real rho4, 
		real rho5){
					matrix[Nx, Ny] K1;
					matrix[Nx, Ny] K2;
					matrix[Nx, Ny] K3;
					matrix[Nx, Ny] Sigma;
	
					//specifying random Gaussian process that governs covariance matrix
					for(i in 1:Nx){
						for (j in 1:Ny){
							K1[i,j] = alpha1*exp(-square(x[i]-y[j])/2/square(rho1));
						}
					}
					
					//specifying random Gaussian process incorporates annual seaonsality (the 1 is a year)
					for(i in 1:Nx){
						for(j in 1:Ny){
							K2[i, j] = alpha2*exp(-2*square(sin(pi()*fabs(x[i]-y[j])*1))/square(rho2))*
							exp(-square(x[i]-y[j])/2/square(rho3));
						}
					}
					
					//specifying random Gaussian process incorporates annual seaonsality (the 1 is a year)
					for(i in 1:Nx){
						for(j in 1:Ny){
							K3[i, j] = alpha3*exp(-square(x[i]-y[j])/2/square(rho4))*
							exp(-square(x[i]-y[j])/2/square(rho5));
						}
					}
						
					Sigma = K1+K2+K3;
					return Sigma;
				}
	//function for posterior calculations
	vector post_pred_rng(
		real a1,
		real a2,
		real a3,
		real r1, 
		real r2,
		real r3,
		real r4, 
		real r5,
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
	Ko = main_GP(No, xo, No, xo, a1, a2, a3, r1, r2, r3, r4, r5);
	for(n in 1:No) Ko[n,n] += sn;
		
	//--------------------------------------------------------------------
	//kernel for predicted data
	Kp = main_GP(Np, xp, Np, xp, a1, a2, a3, r1, r2, r3, r4, r5);
	for(n in 1:Np) Kp[n,n] += sn;
		
	//--------------------------------------------------------------------
	//kernel for observed and predicted cross 
	Kop = main_GP(No, xo, Np, xp, a1, a2, a3, r1, r2, r3, r4, r5);
	
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
	vector[N1] y_log;
	for(n in 1:N1){
	  mu[n] = 0;
	  y_log[n] = log(Y[n]);
	} 
}

parameters {
	real<lower=0> a1;
	real<lower=0> a2;
	real<lower=0> a3;
	real<lower=15> r1;		//Set after some preliminary modeling
	real<lower=0> r2;
	real<lower=0> r3;
	real<lower=0> r4;
	real<lower=0> r5;
	real<lower=0> sigma_sq;
}

model{ 
	matrix[N1,N1] Sigma;
	matrix[N1,N1] L_S;
	
	//using GP function from above 
	Sigma = main_GP(N1, X, N1, X, a1, a2, a3, r1, r2, r3, r4, r5);
	for(n in 1:N1) Sigma[n,n] += sigma_sq;
	
	L_S = cholesky_decompose(Sigma);
	y_log ~ multi_normal_cholesky(mu, L_S);
	
	//priors for parameters
	a1 ~ normal(8.76, 2.65);	//Taken from the third model
	a2 ~ normal(0.05, 0.07);  //Taken from the third model
	a3 ~ student_t(3,0,1);
	r1 ~ normal(24.63, 4.74);	//Taken from the third model
	r2 ~ normal(0.71, 0.13);  //Taken from the third model
	r3 ~ normal(20.92, 9.81);	//Taken from the third model
	r4 ~ student_t(3,0,1);
	r5 ~ student_t(3,0,1);	
	sigma_sq ~ normal(0,1);
}

generated quantities {
	vector[N2] Ypred = post_pred_rng(a1, a2, a3, r1, r2, r3, r4, r5, sigma_sq, N1, X, N2, Xp, y_log);
}