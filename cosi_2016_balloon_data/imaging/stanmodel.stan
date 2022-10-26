data {
    
    int<lower=0> N; // number of response matrix points times number of hours in data set
    int<lower=0> Nh; // number of hours of observation
    int<lower=0> y[Nh*N]; // y observations  
    int<lower=1> Nsky; // number of sky models

    real bg_model[Nh,N]; // BG model
    real conv_sky[Nsky,Nh,N]; // SKY model(s)
    real acceleration_factor_limit;
  
    // background re-normalisation times
    int<lower=0> Ncuts;
    int bg_cuts[Nh];
    int bg_idx_arr[Nh];
  
    //priors on the fitted parameters
    real mu_flux[Nsky];
    real sigma_flux[Nsky];
    real mu_Abg[Ncuts];
    real sigma_Abg[Ncuts];
  
}


//transformed data {
//    // data
//    int data_values[Nh*N];
//  
//    for (nh in 1:Nh) {
//        for (nn in 1:N) {
//            data_values[N*(nh-1)+nn] = y[nh,nn];
//        }
//    }
//
//}


parameters {
    
    real<lower=1e-8, upper=acceleration_factor_limit> flux[Nsky]; // 511 keV line flux of conved sky model
    real<lower=1e-8> Abg[Ncuts]; // background model amplitude(s)
    
}


transformed parameters {

    // model
    real model_values[Nh*N];
    
    for (nh in 1:Nh) {
        for (nn in 1:N) {
            model_values[N*(nh-1)+nn] = Abg[bg_idx_arr[nh]] * bg_model[nh,nn];
            for (ns in 1:Nsky) {
                model_values[N*(nh-1)+nn] += flux[ns] * conv_sky[ns,nh,nn];
            }
        }
    }
    
}



model {

    // normal priors
    flux ~ normal(mu_flux,sigma_flux);
    Abg ~ normal(mu_Abg,sigma_Abg);

    //print(flux);
    //print(Abg);
    //print(model_values);

    // likelihood
    y ~ poisson(model_values);

}


// generated quantities {

//   vector[Nh*N] ppc;
  
//   vector[Nh] model_tot = rep_vector(0,Nh);
//   vector[Nh] model_bg = rep_vector(0,Nh);
//   matrix[Nsky,Nh] model_sky = rep_matrix(0,Nsky,Nh);

//   // create posterior samples for PPC
//   // and
//   // generate the posterior of the
//   for (nh in 1:Nh) {
//       for (nn in 1:N) {
//          ppc[N*(nh-1)+nn] = poisson_rng(model_values[N*(nh-1)+nn]);
//          // fitted model, summed over phi/psi/chi-dimension (only hours left)
//          // sky
//          for (ns in 1:Nsky) {
//              model_sky[ns,nh] += flux[ns] * conv_sky[ns,nh,nn];
//          }
//          // bg
//          model_bg[nh] += Abg[bg_idx_arr[nh]] * bg_model[nh,nn];
//          // total
//          for (ns in 1:Nsky) {
//              model_tot[nh] += flux[ns] * conv_sky[ns,nh,nn];
//          }
//          model_tot[nh] += Abg[bg_idx_arr[nh]] * bg_model[nh,nn];
//       }
//   }

//}