//#include <math.h>
#include <stdio.h>
#include <RcppEigen.h>
using namespace Rcpp;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd pstdnorm(const Eigen::VectorXd v){   ///pnorm0
  Eigen::VectorXd res = Eigen::VectorXd::Zero(v.size());
  for (int i = 0L; i < v.size(); ++i) {
    res(i) = 0.5 * (1L - erf( -v(i) * M_SQRT1_2));
  }
  return(res);
}


// [[Rcpp::depends(RcppEigen)]]
Eigen::VectorXd dpnorm0_dexpsigma(const Eigen::VectorXd v, const double sigma){   ///pnorm0_dsigma
  double sigma2 = std::exp(sigma);
  Eigen::VectorXd res(v.size()),res2(v.size());
  sigma2 = sigma2 * sigma2;
  res = - M_1_SQRT_2PI * v/sigma2;
  res2 = - 0.5 * v.cwiseProduct(v) / sigma2;
  return(res.cwiseProduct(res2.array().exp().matrix()));
}


// [[Rcpp::depends(RcppEigen)]]
Eigen::VectorXd pnorm0(const Eigen::VectorXd v, const double sigma){   ///pnorm0
  Eigen::VectorXd res(v.size());
  double sigma2 = std::exp(sigma);
  for (int i = 0L; i < v.size(); ++i) {
    res(i) = 0.5 * (1L - erf( -v(i) * M_SQRT1_2 / sigma2));
  }
  return(res);
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd dstdnorm(const Eigen::VectorXd v){   ///dnorm
  Eigen::VectorXd v2 = - v.cwiseProduct(v) * 0.5;
  //Eigen::VectorXd res = 0.5 * M_SQRT1_2 * M_2_SQRTPI * v2.array().exp().matrix();
  Eigen::VectorXd res = M_1_SQRT_2PI * v2.array().exp().matrix();
  return(res);
}


// [[Rcpp::depends(RcppEigen)]]
Eigen::VectorXd dnorm0(const Eigen::VectorXd v, const double sigma){   ///dnorm0
  double s1 = 1 / std::exp(sigma);
  Eigen::VectorXd v2 = v * s1;
  v2 = - v2.cwiseProduct(v2) * 0.5;
  //Eigen::VectorXd res = 0.5 * M_SQRT1_2 * M_2_SQRTPI * v2.array().exp().matrix() * s1;
  Eigen::VectorXd res = M_1_SQRT_2PI * v2.array().exp().matrix() * s1;
  return(res);
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd dstdlogis(const Eigen::VectorXd v) {
  Eigen::VectorXd res = v.array().exp().matrix();
  Eigen::VectorXd res1 = res + res.Ones(res.size());
  return(res.cwiseQuotient(res1.cwiseProduct(res1)));
}


// [[Rcpp::depends(RcppEigen)]]
Eigen::VectorXd dlogis0(const Eigen::VectorXd v, const double s) {
  double s2 = std::exp(s);
  Eigen::VectorXd res = (v / s2).array().exp().matrix();
  Eigen::VectorXd res1 = res + res.Ones(res.size());
  return(res.cwiseQuotient(res1.cwiseProduct(res1)) / s2);
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd pstdlogis(const Eigen::VectorXd v) {
  Eigen::VectorXd res = (-v).array().exp().matrix();
  return((res + res.Ones(res.size())).cwiseInverse());
}


// [[Rcpp::depends(RcppEigen)]]
Eigen::VectorXd plogis0(const Eigen::VectorXd v, const double s) {
  Eigen::VectorXd res = (-v / std::exp(s)).array().exp().matrix();
  return((res + res.Ones(res.size())).cwiseInverse());
}


// [[Rcpp::depends(RcppEigen)]]
Eigen::VectorXd dplogis0_dexps(const Eigen::VectorXd v, const double s) {
  double s2 = std::exp(s);
  Eigen::VectorXd res = (-v / s2).array().exp().matrix();
  Eigen::VectorXd res1 =  - res.cwiseProduct(v / s2 /s2);
  res = res + res.Ones(res.size());
  res = res.cwiseProduct(res);
  return(res1.cwiseQuotient(res));
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd colpath(const Eigen::MatrixXd m,
                        const Eigen::VectorXi v,  // indexing start from 1 (R)
                        const int offset) {
  Eigen::VectorXd res = Eigen::VectorXd::Zero(v.size());
  for (int i = 0L; i < v.size(); ++i) {
    res(i) = m(i, v(i) - 1L + offset);
  }
  return(res);
}


// [[Rcpp::depends(RcppEigen)]]
Eigen::VectorXd extract_elements(const Eigen::VectorXi x,
                                 const int offset,
                                 const Eigen::VectorXd v){
  Eigen::VectorXd res = Eigen::VectorXd::Zero(x.size());
  for (int i = 0L; i < x.size(); ++i) {
    res(i) = v( x(i) + offset );
  }
  return(res);
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
SEXP glm2hopit(const Eigen::VectorXd latent_params,
               const Eigen::VectorXd thresh_lambda,
               const Eigen::VectorXd thresh_gamma,
               const int thresh_1_exp){

  int J_1 = thresh_lambda.size();

  Eigen::VectorXd latent_params_g = latent_params;  // already changed to negative

  Eigen::VectorXd thresh_lambda_g = thresh_lambda;
  Eigen::MatrixXd thresh_gamma_g2 = thresh_gamma;
  thresh_gamma_g2.resize(J_1, thresh_gamma.size()/ J_1); // gamma as J_1 x n_thresh_var matrix
  Eigen::MatrixXd thresh_gamma_g3 = thresh_gamma_g2;

  if (thresh_1_exp == 1) {
    thresh_lambda_g(0) = 0; //cannnot take log from negative number
  }
  for (int i = 1; i < (J_1); ++i){ // in glm lambda is cumulated
    thresh_lambda_g(i) = std::log(thresh_lambda(i) - thresh_lambda(i-1));
  }
  for (int i = 1; i < (J_1); ++i){   // in glm gamma equvalent is cumulated
    thresh_gamma_g2.row(i) = thresh_gamma_g3.row(i) - thresh_gamma_g3.row(i-1);
  }

  Eigen::VectorXd thresh_gamma_g(Eigen::Map<Eigen::VectorXd>(thresh_gamma_g2.data(), thresh_gamma.size()));
  Eigen::VectorXd coef(latent_params.size()+thresh_lambda.size()+thresh_gamma.size());
  coef << latent_params_g, thresh_lambda_g, thresh_gamma_g;
  return List::create(
    Named("latent_params") = latent_params_g,
    Named("thresh_lambda") = thresh_lambda_g,
    Named("thresh_gamma") = thresh_gamma_g,
    Named("coef") = coef
  );
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
SEXP glm2hopit_nogamma(const Eigen::VectorXd latent_params,
                       const Eigen::VectorXd thresh_lambda,
                       const int thresh_1_exp){

  int J_1 = thresh_lambda.size();

  Eigen::VectorXd latent_params_g = latent_params;  // already changed to negative
  Eigen::VectorXd thresh_lambda_g = thresh_lambda;

  if (thresh_1_exp == 1) {
    thresh_lambda_g(0) = 0; //cannnot take log from negative number
  }
  for (int i = 1; i < (J_1); ++i){ // in glm lambda is cumulated
    thresh_lambda_g(i) = std::log(thresh_lambda(i) - thresh_lambda(i-1));
  }

  Eigen::VectorXd coef(latent_params.size()+thresh_lambda.size());
  coef << latent_params_g, thresh_lambda_g;
  return List::create(
    Named("latent_params") = latent_params_g,
    Named("thresh_lambda") = thresh_lambda_g,
    Named("coef") = coef
  );
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd ind_latent_thresh(const Eigen::MatrixXd thresh_mm,
                                  const Eigen::VectorXd thresh_lambda,
                                  const Eigen::VectorXd thresh_gamma){
  int J_1 = thresh_lambda.size();
  Eigen::MatrixXd res(thresh_mm.rows(), J_1);
  for (int i = 0L; i < thresh_mm.rows(); ++i) {
    res.row(i) = thresh_lambda;
  }

  Eigen::VectorXi pos = Eigen::VectorXi::Zero(J_1);
  for (int i = 0L; i < J_1; ++i) pos[i] = i;

  for (int i = 0L; i < thresh_mm.cols(); ++i) {
    int offse = i * J_1;
    res = res + thresh_mm.col(i) * extract_elements(pos, offse, thresh_gamma).transpose();
  }
  return(res);
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd getThresholds(
    const Eigen::MatrixXd thresh_mm,
    const Eigen::VectorXd thresh_lambda,
    const Eigen::VectorXd thresh_gamma,
    const int thresh_no_cov,
    const double thresh_start,  // zero threshold , typically -Inf
    const int thresh_1_exp) {  // if to exponentiate threshold 1

  int J_1 = thresh_lambda.size();
  long N = thresh_mm.rows();
  Eigen::VectorXd t_thresh_lambda = thresh_lambda.transpose();

  Eigen::MatrixXd Lin_Thresh_mat(N, J_1);
  if (thresh_no_cov == 1L) {
    for (long i = 0L; i < N; ++i) {     // maybe fill by columns to increase speed
      Lin_Thresh_mat.row(i) = thresh_lambda;
    }
  } else {
    Lin_Thresh_mat = ind_latent_thresh(thresh_mm, thresh_lambda, thresh_gamma);
  }

  Eigen::MatrixXd a = Eigen::MatrixXd::Zero(N, J_1 + 2L);
  a.col(J_1 + 1L).fill(R_PosInf);
  a.col(0L).fill(thresh_start);
  if (thresh_1_exp == 0L) a.col(1L) = Lin_Thresh_mat.col(0L); else a.col(1L) = Lin_Thresh_mat.col(0L).array().exp().matrix();

  Lin_Thresh_mat = Lin_Thresh_mat.array().exp().matrix();

  for  (int i = 2L; i <= J_1; ++i)  a.col(i) = a.col(i - 1L) + Lin_Thresh_mat.col(i - 1L);
  return(a);
}


// [[Rcpp::depends(RcppEigen)]]
Eigen::VectorXd subvec(const int sta, const int end, const Eigen::VectorXd v){  // this subseting is slightlyfaster than next one, why???
  Eigen::VectorXd res = Eigen::VectorXd::Zero(end-sta + 1L);
  int j = 0;
  for  (int i = sta; i <= end; ++i) {
    res(j) = v(i);
    j++;
  }
  return(res);
}

// [[Rcpp::depends(RcppEigen)]]
Eigen::VectorXd subvec_b(const int sta, const int end, const Eigen::VectorXd v){
  return(v.segment(sta,end-sta+1));
}

// [[Rcpp::depends(RcppEigen)]]
Eigen::VectorXd exchange(const Eigen::VectorXd x, const double from, const double to){
  Eigen::VectorXd res(x.size());
  for  (long i = 0; i < x.size(); ++i) {
    if (x(i) == from) {
      res(i) = to;
    } else res(i) = x(i);
  }
  return(res);
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double LLFunc(const Eigen::Map<Eigen::VectorXd> parameters,
              const Eigen::VectorXi yi,
              const Eigen::MatrixXd latent_mm,
              const Eigen::MatrixXd thresh_mm,
              const Eigen::VectorXi parcount,
              const int hasdisp, // 0-no, 1 - yes
              const int link, /// 0-probit , 1 - logit
              const int thresh_no_cov, /// 0-NO , 1-Yes
              const int negative, ///1 -yes 0 no
              const int use_weights,
              const Eigen::VectorXd weights,
              const double thresh_start,
              const int thresh_1_exp,
              const double out_val){ /// e.g. -inf

  Eigen::VectorXd latent_par = subvec(0, parcount(0) - 1L, parameters);
  int p01 = parcount(0) + parcount(1) - 1L;
  Eigen::VectorXd thresh_lambda = subvec(parcount(0), p01, parameters);
  Eigen::VectorXd thresh_gamma = subvec(p01 + 1L, p01 + parcount(2), parameters);
  double disp = 0;
  if (hasdisp == 1) disp = parameters(parameters.size()-1);

  Eigen::MatrixXd a;
  Eigen::MatrixXd b;
  a = getThresholds(thresh_mm, thresh_lambda, thresh_gamma, thresh_no_cov, thresh_start, thresh_1_exp);
  b = latent_mm * latent_par;

  Eigen::VectorXd LO = colpath(a, yi, 0L) - b;
  Eigen::VectorXd HI = colpath(a, yi, 1L) - b;
  LO = exchange(LO, R_NegInf, -20);
  HI = exchange(HI, R_PosInf, 20);

  Eigen::VectorXd P;
  if (hasdisp == 1){
    if (link == 0)  P = pnorm0(HI, disp) - pnorm0(LO, disp); else P = plogis0(HI, disp) - plogis0(LO, disp);
  } else {
    if (link == 0)  P = pstdnorm(HI) - pstdnorm(LO); else P = pstdlogis(HI) - pstdlogis(LO);
  }

  int d;
  if (negative == 1L) d = -1; else d = 1;

  double res;

  if (P.minCoeff() > 0L) {
    P = d * P.array().log().matrix();
    if (use_weights == 1) P = P.cwiseProduct(weights);
    res = P.array().sum();
  } else {
    if (out_val == R_NegInf) res = out_val * d; else {
      double out_val2 = std::exp(out_val);
      for (int i=0; i<P.size(); ++i) if (P(i)<=0) P(i) = out_val2;
      P = d * P.array().log().matrix();
      res = P.array().sum();
    }
  }
  return(res);
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd LLFuncIndv(const Eigen::Map<Eigen::VectorXd> parameters,
                           const Eigen::VectorXi yi,
                           const Eigen::MatrixXd latent_mm,
                           const Eigen::MatrixXd thresh_mm,
                           const Eigen::VectorXi parcount,
                           const int hasdisp, // 0-no, 1 - yes
                           const int link, /// 0-probit , 1 - logit
                           const int thresh_no_cov, /// 0-NO , 1-Yes
                           const int negative, ///1 -yes 0 no
                           const int use_weights,
                           const double thresh_start,
                           const int thresh_1_exp,
                           const Eigen::VectorXd weights){ /// 1-yes, 0-no

  Eigen::VectorXd latent_par = subvec(0, parcount(0) - 1L, parameters);
  int p01 = parcount(0) + parcount(1) - 1L;
  Eigen::VectorXd thresh_lambda = subvec(parcount(0), p01, parameters);
  Eigen::VectorXd thresh_gamma = subvec(p01 + 1L, p01 + parcount(2), parameters);
  double disp = 0;
  if (hasdisp == 1) disp = parameters(parameters.size()-1);
  Eigen::MatrixXd a;
  Eigen::MatrixXd b;
  a = getThresholds(thresh_mm, thresh_lambda, thresh_gamma, thresh_no_cov, thresh_start, thresh_1_exp);
  b = latent_mm * latent_par;

  Eigen::VectorXd LO = colpath(a, yi, 0L) - b;
  Eigen::VectorXd HI = colpath(a, yi, 1L) - b;
  LO = exchange(LO, R_NegInf, -20);
  HI = exchange(HI, R_PosInf, 20);

  Eigen::VectorXd P;

  if (hasdisp == 1){
    if (link == 0)  P = pnorm0(HI, disp) - pnorm0(LO, disp); else P = plogis0(HI, disp) - plogis0(LO, disp);
  } else {
    if (link == 0)  P = pstdnorm(HI) - pstdnorm(LO); else P = pstdlogis(HI) - pstdlogis(LO);
  }

  int d;
  if (negative == 1L) d = -1;  else d = 1;
  P = d * P.array().log().matrix();
  if (use_weights == 1) P = P.cwiseProduct(weights);
  return(P);
}


// [[Rcpp::depends(RcppEigen)]]
Eigen::MatrixXd rep_col_c(const Eigen::MatrixXd mat, int times) {
  Eigen::MatrixXd res(mat.rows(), mat.cols() * times);
  int k = 0;
  for (int i=0; i< times; ++i){
    for (int j=0; j< mat.cols(); ++j){
      res.col(k)=mat.col(j);
      k++;
    }
  }
  return(res);
}


// [[Rcpp::depends(RcppEigen)]]
Eigen::VectorXd colSums_c(const Eigen::MatrixXd mat) {
  Eigen::VectorXd res(mat.cols());
  for (int i=0; i< mat.cols(); ++i){
    res(i) = mat.col(i).array().sum();
  }
  return (res);
}


// [[Rcpp::depends(RcppEigen)]]
Eigen::VectorXd colSums_weighted(const Eigen::MatrixXd mat, const Eigen::VectorXd weights) {
  Eigen::VectorXd res(mat.cols());
  for (int i=0; i< mat.cols(); ++i){
    res(i) = (mat.col(i).cwiseProduct(weights)).array().sum();
  }
  return (res);
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd LLGradFunc(const Eigen::Map<Eigen::VectorXd> parameters,
                           const Eigen::VectorXi yi,
                           const Eigen::MatrixXd YYY1,
                           const Eigen::MatrixXd YYY2,
                           const Eigen::MatrixXd YYY3,
                           const Eigen::MatrixXd YYY4,
                           const Eigen::MatrixXd latent_mm,
                           const Eigen::MatrixXd thresh_mm,
                           const Eigen::MatrixXd thresh_extd,
                           const Eigen::VectorXi parcount,
                           const int hasdisp, // 0-no, 1 - yes
                           const int link, /// 0-probit , 1 - logit
                           const int thresh_no_cov, /// 0-NO , 1-Yes
                           const int negative, ///1 -yes 0 no
                           const int use_weights, ///1 -yes 0 no
                           const double thresh_start,
                           const int thresh_1_exp,
                           const Eigen::VectorXd weights){ /// 1-yes, 0-no

  long N = yi.size();

  Eigen::VectorXd latent_par = subvec(0, parcount(0) - 1L, parameters);
  int p01 = parcount(0) + parcount(1) - 1L;
  Eigen::VectorXd thresh_lambda = subvec(parcount(0), p01, parameters);
  Eigen::VectorXd thresh_gamma = subvec(p01 + 1L, p01 + parcount(2), parameters);
  double disp = 0;
  if (hasdisp == 1) disp = parameters(parameters.size()-1);

  Eigen::MatrixXd a;
  Eigen::MatrixXd b;
  a = getThresholds(thresh_mm, thresh_lambda, thresh_gamma, thresh_no_cov, thresh_start, thresh_1_exp);
  b = latent_mm * latent_par;

  Eigen::VectorXd A2 = colpath(a, yi, 0L) - b;
  Eigen::VectorXd A1 = colpath(a, yi, 1L) - b;
  A2 = exchange(A2, R_NegInf, -20);
  A1 = exchange(A1, R_PosInf, 20);

  Eigen::VectorXd D, P, P1, P2, D1, D2, sig1, sig2;
  if (hasdisp == 1){
    if (link == 0) {
      P1 = pnorm0(A1,disp);
      P2 = pnorm0(A2,disp);
      D1 = dnorm0(A1,disp);
      D2 = dnorm0(A2,disp);
      sig1 = dpnorm0_dexpsigma(A1, disp);
      sig2 = dpnorm0_dexpsigma(A2, disp);
    } else {
      P1 = plogis0(A1,disp);
      P2 = plogis0(A2,disp);
      D1 = dlogis0(A1,disp);
      D2 = dlogis0(A2,disp);
      sig1 = dplogis0_dexps(A1, disp);
      sig2 = dplogis0_dexps(A2, disp);
    }
  } else {
    if (link == 0) {
      P1 = pstdnorm(A1);//pstdnorm(A1);
      P2 = pstdnorm(A2);//pstdnorm(A2);
      D1 = dstdnorm(A1);
      D2 = dstdnorm(A2);
    } else {
      P1 = pstdlogis(A1);
      P2 = pstdlogis(A2);
      D1 = dstdlogis(A1);
      D2 = dstdlogis(A2);
    }
  }

  P = P1 - P2;
  D = D1 - D2;

  Eigen::VectorXd dlnLL_dX = P.array().cwiseInverse().matrix();

  Eigen::VectorXd dF_dexps;
  if (hasdisp == 1) dF_dexps = (sig1 - sig2).cwiseProduct(dlnLL_dX);

  Eigen::MatrixXd dlnLL_dbeta(D.size(), latent_par.size());
  Eigen::MatrixXd dd = - D.cwiseProduct(dlnLL_dX);
  for (int i = 0L; i < latent_par.size(); ++i) {
    dlnLL_dbeta.col(i) = dd;
  }
  dlnLL_dbeta = dlnLL_dbeta.cwiseProduct(latent_mm);

  Eigen::MatrixXd da;
  Eigen::MatrixXd dlnLL_Lambda(D1.size(), a.cols() - 2);
  Eigen::MatrixXd dlnLL_Gamma;
  Eigen::VectorXd res;

  Eigen::MatrixXd D1_(D1.rows(), YYY1.cols()), D2_(D1.rows(), YYY1.cols()), D3_(D1.rows(), YYY1.cols());
  da = Eigen::MatrixXd::Ones(N, YYY1.cols());

  for (int i = 1; i < (a.cols() - 2); ++i){
    da.col(i) = a.col(i + 1) - a.col(i);
  }

  for (int i = 0; i < da.cols(); ++i){
    D1_.col(i) = D1;
    D2_.col(i) = D2;
    D3_.col(i) = dlnLL_dX;
  }
  dlnLL_Lambda = (D1_.cwiseProduct(da.cwiseProduct(YYY1)) - D2_.cwiseProduct(da.cwiseProduct(YYY2))).cwiseProduct(D3_);

  if (use_weights==0){
    if (thresh_no_cov == 1L) {
      if (hasdisp == 1) {
        res.resize(dlnLL_dbeta.cols()+dlnLL_Lambda.cols()+1);
        res << colSums_c(dlnLL_dbeta), colSums_c(dlnLL_Lambda), colSums_c(dF_dexps);
      } else {
        res.resize(dlnLL_dbeta.cols()+dlnLL_Lambda.cols());
        res << colSums_c(dlnLL_dbeta), colSums_c(dlnLL_Lambda);
      }
    } else {
      dlnLL_Gamma = rep_col_c(dlnLL_Lambda, thresh_mm.cols()).cwiseProduct(thresh_extd);
      if (hasdisp == 1) {
        res.resize(dlnLL_dbeta.cols()+dlnLL_Lambda.cols()+dlnLL_Gamma.cols()+1);
        res << colSums_c(dlnLL_dbeta), colSums_c(dlnLL_Lambda), colSums_c(dlnLL_Gamma), colSums_c(dF_dexps);
      } else {
        res.resize(dlnLL_dbeta.cols()+dlnLL_Lambda.cols()+dlnLL_Gamma.cols());
        res << colSums_c(dlnLL_dbeta), colSums_c(dlnLL_Lambda), colSums_c(dlnLL_Gamma);
      }
    }
  } else {
    if (thresh_no_cov == 1L) {
      if (hasdisp == 1) {
        res.resize(dlnLL_dbeta.cols()+dlnLL_Lambda.cols() +1);
        res << colSums_weighted(dlnLL_dbeta, weights), colSums_weighted(dlnLL_Lambda, weights), colSums_weighted(dF_dexps, weights); // weigth sums by weights
      } else {
        res.resize(dlnLL_dbeta.cols()+dlnLL_Lambda.cols());
        res << colSums_weighted(dlnLL_dbeta, weights), colSums_weighted(dlnLL_Lambda, weights); // weigth sums by weights
      }
    } else {
      dlnLL_Gamma = rep_col_c(dlnLL_Lambda, thresh_mm.cols()).cwiseProduct(thresh_extd);
      if (hasdisp == 1) {
        res.resize(dlnLL_dbeta.cols()+dlnLL_Lambda.cols()+dlnLL_Gamma.cols()+1);
        res << colSums_weighted(dlnLL_dbeta, weights), colSums_weighted(dlnLL_Lambda, weights),
               colSums_weighted(dlnLL_Gamma, weights), colSums_weighted(dF_dexps, weights);
      } else {
        res.resize(dlnLL_dbeta.cols()+dlnLL_Lambda.cols()+dlnLL_Gamma.cols());
        res << colSums_weighted(dlnLL_dbeta, weights), colSums_weighted(dlnLL_Lambda, weights), colSums_weighted(dlnLL_Gamma, weights);
      }
    }
  }
  if (negative == 1) res = - res;
  return(res);
}


// [[Rcpp::depends(RcppEigen)]]
Eigen::MatrixXd weight_rows(const Eigen::MatrixXd m, const Eigen::VectorXd weights){
  Eigen::MatrixXd res(m);
  for (int i=0; i< m.cols(); ++i){
    res.col(i) = m.col(i).cwiseProduct(weights);
  }
  return (res);
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd LLGradFuncIndv(const Eigen::Map<Eigen::VectorXd> parameters,
                               const Eigen::VectorXi yi,
                               const Eigen::MatrixXd YYY1,
                               const Eigen::MatrixXd YYY2,
                               const Eigen::MatrixXd YYY3,
                               const Eigen::MatrixXd YYY4,
                               const Eigen::MatrixXd latent_mm,
                               const Eigen::MatrixXd thresh_mm,
                               const Eigen::MatrixXd thresh_extd,
                               const Eigen::VectorXi parcount,
                               const int hasdisp, // 0-no, 1 - yes
                               const int link, /// 0-probit , 1 - logit
                               const int thresh_no_cov, /// 0-NO , 1-Yes
                               const int negative, ///1 -yes 0 no
                               const int use_weights, ///1 -yes 0 no
                               const double thresh_start,
                               const int thresh_1_exp,
                               const Eigen::VectorXd weights){ /// 1-yes, 0-no

  long N = yi.size();

  Eigen::VectorXd latent_par = subvec(0, parcount(0) - 1L, parameters);
  int p01 = parcount(0) + parcount(1) - 1L;
  Eigen::VectorXd thresh_lambda = subvec(parcount(0), p01, parameters);
  Eigen::VectorXd thresh_gamma = subvec(p01 + 1L, p01 + parcount(2), parameters);
  double disp = 0;
  if (hasdisp == 1) disp = parameters(parameters.size()-1);

  Eigen::MatrixXd a;
  Eigen::MatrixXd b;
  a = getThresholds(thresh_mm, thresh_lambda, thresh_gamma, thresh_no_cov, thresh_start, thresh_1_exp);
  b = latent_mm * latent_par;

  Eigen::VectorXd A2 = colpath(a, yi, 0L) - b;
  Eigen::VectorXd A1 = colpath(a, yi, 1L) - b;
  A2 = exchange(A2, R_NegInf, -20);
  A1 = exchange(A1, R_PosInf, 20);

  Eigen::VectorXd D, P, P1, P2, D1, D2, sig1, sig2;
  if (hasdisp == 1){
    if (link == 0) {
      P1 = pnorm0(A1,disp);
      P2 = pnorm0(A2,disp);
      D1 = dnorm0(A1,disp);
      D2 = dnorm0(A2,disp);
      sig1 = dpnorm0_dexpsigma(A1, disp);
      sig2 = dpnorm0_dexpsigma(A2, disp);
    } else {
      P1 = plogis0(A1,disp);
      P2 = plogis0(A2,disp);
      D1 = dlogis0(A1,disp);
      D2 = dlogis0(A2,disp);
      sig1 = dplogis0_dexps(A1, disp);
      sig2 = dplogis0_dexps(A2, disp);
    }
  } else {
    if (link == 0) {
      P1 = pstdnorm(A1);//pstdnorm(A1);
      P2 = pstdnorm(A2);//pstdnorm(A2);
      D1 = dstdnorm(A1);
      D2 = dstdnorm(A2);
    } else {
      P1 = pstdlogis(A1);
      P2 = pstdlogis(A2);
      D1 = dstdlogis(A1);
      D2 = dstdlogis(A2);
    }
  }

  P = P1 - P2;
  D = D1 - D2;

  Eigen::VectorXd dlnLL_dX = P.array().cwiseInverse().matrix();

  Eigen::VectorXd dF_dexps;
  if (hasdisp == 1) dF_dexps = (sig1 - sig2).cwiseProduct(dlnLL_dX);

  Eigen::MatrixXd dlnLL_dbeta = Eigen::MatrixXd::Zero(D.size(), latent_par.size());
  Eigen::MatrixXd dd = - D.cwiseProduct(dlnLL_dX);
  for (int i = 0L; i < latent_par.size(); ++i) {
    dlnLL_dbeta.col(i) = dd;
  }

  dlnLL_dbeta = dlnLL_dbeta.cwiseProduct(latent_mm);

  Eigen::MatrixXd da;
  Eigen::MatrixXd dlnLL_Lambda(D1.size(), a.cols() - 2);
  Eigen::MatrixXd dlnLL_Gamma;

  Eigen::MatrixXd D1_(D1.rows(), YYY1.cols()), D2_(D1.rows(), YYY1.cols()), D3_(D1.rows(), YYY1.cols());
  da = Eigen::MatrixXd::Ones(N, YYY1.cols());

  for (int i = 1; i < (a.cols() - 2); ++i){
    da.col(i) = a.col(i + 1) - a.col(i);
  }

  for (int i = 0; i < da.cols(); ++i){
    D1_.col(i) = D1;
    D2_.col(i) = D2;
    D3_.col(i) = dlnLL_dX;
  }
  dlnLL_Lambda = (D1_.cwiseProduct(da.cwiseProduct(YYY1)) - D2_.cwiseProduct(da.cwiseProduct(YYY2))).cwiseProduct(D3_);

  int k;
  if (thresh_no_cov == 1L) k = 0; else k = parcount(2);
  if (hasdisp == 1) k++;
  Eigen::MatrixXd res(dlnLL_dbeta.rows(), dlnLL_dbeta.cols() + dlnLL_Lambda.cols() + k);

  if (use_weights==0){
    if (thresh_no_cov == 1L) {
      if (hasdisp == 1) {
        res << dlnLL_dbeta, dlnLL_Lambda, dF_dexps;
      } else {
        res << dlnLL_dbeta, dlnLL_Lambda;
      }
    } else {
      Eigen::MatrixXd dlnLL_Gamma = rep_col_c(dlnLL_Lambda, thresh_mm.cols()).cwiseProduct(thresh_extd);
      if (hasdisp == 1) {
        res << dlnLL_dbeta, dlnLL_Lambda, dlnLL_Gamma, dF_dexps;
      } else {
        res << dlnLL_dbeta, dlnLL_Lambda, dlnLL_Gamma;
      }
    }
  } else {
    if (thresh_no_cov == 1L) {
      if (hasdisp == 1) {
        res << weight_rows(dlnLL_dbeta, weights), weight_rows(dlnLL_Lambda, weights), dF_dexps;
      } else {
        res << weight_rows(dlnLL_dbeta, weights), weight_rows(dlnLL_Lambda, weights); // weigth sums by weights
      }
    } else {
      Eigen::MatrixXd dlnLL_Gamma = rep_col_c(dlnLL_Lambda, thresh_mm.cols()).cwiseProduct(thresh_extd);
      if (hasdisp == 1) {
        res << weight_rows(dlnLL_dbeta, weights), weight_rows(dlnLL_Lambda, weights), weight_rows(dlnLL_Gamma, weights), dF_dexps;
      } else {
        res << weight_rows(dlnLL_dbeta, weights), weight_rows(dlnLL_Lambda, weights), weight_rows(dlnLL_Gamma, weights);
      }
    }
  }

  if (negative == 1) res = - res;
  return(res);
}
