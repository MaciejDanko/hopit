//#include <math.h>
#include <RcppEigen.h>
using namespace Rcpp;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// [[Rcpp::depends(RcppEigen)]]
Eigen::VectorXd pstdnorm(const Eigen::VectorXd v){   ///pnorm
  Eigen::VectorXd res = Eigen::VectorXd::Zero(v.size());
  for (int i = 0L; i < v.size(); ++i) {
    res(i) = 0.5 * (1L - erf( -v(i) * M_SQRT1_2));
  }
  return(res);
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd dstdnorm(const Eigen::VectorXd v){   ///dnorm
  Eigen::VectorXd v2 = - v.cwiseProduct(v) * 0.5;
  Eigen::VectorXd res = 0.5 * M_SQRT1_2 * M_2_SQRTPI * v2.array().exp().matrix();
  return(res);
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd pstdlogit(const Eigen::VectorXd v) {
  Eigen::VectorXd res = v.array().exp().matrix();
  Eigen::VectorXd res1 = res + res.Ones(res.size());
  return(res.cwiseQuotient(res1));
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd dstdlogit(const Eigen::VectorXd v) {
  Eigen::VectorXd res = (-v).array().exp().matrix();
  Eigen::VectorXd res1 = res + res.Ones(res.size());
  res1 = res1.cwiseProduct(res1);
  return(res.cwiseQuotient(res1));
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

// [[Rcpp::export]]
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
SEXP vglm2gotm_jurges_exp(const Eigen::VectorXd reg_params,
                          const Eigen::VectorXd thresh_lambda,
                          const Eigen::VectorXd thresh_gamma,
                          const int thresh_fun){
  int J_1 = thresh_lambda.size();
  Eigen::VectorXd reg_params_g = - reg_params;
  Eigen::VectorXd thresh_lambda_g = thresh_lambda;
  Eigen::MatrixXd thresh_gamma_g2 = thresh_gamma;
  thresh_gamma_g2.resize(J_1, thresh_gamma.size()/ J_1);
  Eigen::MatrixXd thresh_gamma_g3 = thresh_gamma_g2;
  if (thresh_fun==0) {
    for (int i = 1; i < (J_1); ++i){
      thresh_lambda_g(i) = std::log(thresh_lambda(i) - thresh_lambda(i-1));
    }
    for (int i = 1; i < (J_1); ++i){
      thresh_gamma_g2.row(i) = thresh_gamma_g3.row(i) - thresh_gamma_g3.row(i-1);
    }
  } else {
    for (int i = 1; i < (J_1); ++i){
      thresh_lambda_g(i) = std::log(thresh_lambda(i) - thresh_lambda(i-1));
    }
    for (int i = 1; i < (J_1); ++i){
      thresh_gamma_g2.row(i) = (thresh_gamma_g3.row(i) - thresh_gamma_g3.row(i-1)).array().exp().matrix();
    }
  }

  Eigen::VectorXd thresh_gamma_g(Eigen::Map<Eigen::VectorXd>(thresh_gamma_g2.data(), thresh_gamma.size()));
  Eigen::VectorXd coef(reg_params.size()+thresh_lambda.size()+thresh_gamma.size());
  coef << reg_params_g, thresh_lambda_g, thresh_gamma_g;
  return List::create(
    Named("reg_params") = reg_params_g,
    Named("thresh_lambda") = thresh_lambda_g,
    Named("thresh_gamma") = thresh_gamma_g,
    Named("coef") = coef
  );
}


// Eigen::SparseMatrix<double> diff(Eigen::SparseMatrix<double> E) {
//   Eigen::SparseMatrix<double> E1 = E.block(0, 0, E.rows()-1, E.cols());
//   Eigen::SparseMatrix<double> E2 = E.block(1, 0, E.rows()-1, E.cols());
//   return E2 - E1;
// }

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd ind_reg_thresh(const Eigen::MatrixXd thresh_mm,
                               const Eigen::VectorXd thresh_lambda,
                               const Eigen::VectorXd thresh_gamma){
  int J_1 = thresh_lambda.size();
  Eigen::MatrixXd res = Eigen::MatrixXd::Zero(thresh_mm.rows(), J_1);
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
    const int thresh_no_cov, /// 0-NO , 1-Yes
    const int thresh_method, /// 0-Jurges, 1-Hopit
    const int thresh_func, /// 0-exp, 1-lin
    const int use_alpha, /// 0-auto
    const double alpha_0 ){ ///

  int J_1 = thresh_lambda.size();
  long N = thresh_mm.rows();
  Eigen::VectorXd t_thresh_lambda = thresh_lambda.transpose();

  Eigen::MatrixXd Lin_Thresh_mat;
  if (thresh_no_cov == 1L) {
    Lin_Thresh_mat = Eigen::MatrixXd::Zero(N, J_1);
    for (long i = 0L; i < N; ++i) {     // maybe fill by columns to increase speed
      Lin_Thresh_mat.row(i) = thresh_lambda;
    }
  } else {
    Lin_Thresh_mat = ind_reg_thresh(thresh_mm, thresh_lambda, thresh_gamma);
  }

  Eigen::MatrixXd a = Eigen::MatrixXd::Zero(N, J_1 + 2L);
  a.col(J_1 + 1L).fill(R_PosInf); /// last column
  if (thresh_method == 0L){
    if (use_alpha == 1L) {
      a.col(0L).fill(alpha_0);
    } else {
      a.col(0L).fill(R_NegInf);
    }
    a.col(1L) = Lin_Thresh_mat.col(0L);
  } else {
    if (use_alpha == 1L) {
      a.col(0L).fill(alpha_0);
    } else {
      a.col(0L).fill(0L);
    }
    if (thresh_func == 1) {
      a.col(1L) = a.col(0L) + Lin_Thresh_mat.col(0L).array().abs().matrix();  // must be abs here
    } else {
      a.col(1L) = a.col(0L) + Lin_Thresh_mat.col(0L).array().exp().matrix(); // use exp1 if does not work
    }
  }

  if (thresh_func == 0) {
    Lin_Thresh_mat = Lin_Thresh_mat.array().exp().matrix();
  }

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
  Eigen::VectorXd res = Eigen::VectorXd::Zero(x.size());
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
              const Eigen::MatrixXd reg_mm,
              const Eigen::MatrixXd thresh_mm,
              const Eigen::VectorXi parcount,
              const int link, /// 0-probit , 1 - logit
              const int thresh_no_cov, /// 0-NO , 1-Yes
              const int thresh_method, /// 0-Jurges, 1-Hopit
              const int thresh_func, /// 0-exp, 1-lin
              const int use_alpha, /// 0-auto
              const double alpha_0,
              const int negative, ///1 -yes 0 no
              const int use_weights,
              const Eigen::VectorXd weights,
              const double out_val){ /// e.g. -inf

  Eigen::VectorXd reg_par = subvec(0, parcount(0) - 1L, parameters);
  int p01 = parcount(0) + parcount(1) - 1L;
  Eigen::VectorXd thresh_lambda = subvec(parcount(0), p01, parameters);
  Eigen::VectorXd thresh_gamma = subvec(p01 + 1L, p01 + parcount(2), parameters);
  Eigen::MatrixXd a = getThresholds(thresh_mm, thresh_lambda, thresh_gamma, thresh_no_cov,
                                    thresh_method, thresh_func, use_alpha, alpha_0 );
  Eigen::VectorXd b = reg_mm * reg_par;

  Eigen::VectorXd LO = colpath(a, yi, 0L) - b;
  Eigen::VectorXd HI = colpath(a, yi, 1L) - b;
  LO = exchange(LO, R_NegInf, -20);
  HI = exchange(HI, R_PosInf, 20);

  Eigen::VectorXd P;
  if (link == 0) {
    P = pstdnorm(HI) - pstdnorm(LO);
  } else {
    P = pstdlogit(HI) - pstdlogit(LO);
  }
  int d;
  if (negative == 1L) {
    d = -1;
  }  else d = 1;

  double res;

  if (P.minCoeff() > 0L) {
    P = d * P.array().log().matrix();
    if (use_weights == 1) P = P.cwiseProduct(weights);
    res = P.array().sum();
  } else {
    if (out_val == R_NegInf) res = out_val * d; else {
      double out_val2 = std::exp(out_val);
      for (int i=0; i<P.size(); ++i){
        if (P(i)<=0) P(i) = out_val2;
      }
      P = d * P.array().log().matrix();
    }
  }
  return(res);
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd LLFuncIndv(const Eigen::Map<Eigen::VectorXd> parameters,
                           const Eigen::VectorXi yi,
                           const Eigen::MatrixXd reg_mm,
                           const Eigen::MatrixXd thresh_mm,
                           const Eigen::VectorXi parcount,
                           const int link, /// 0-probit , 1 - logit
                           const int thresh_no_cov, /// 0-NO , 1-Yes
                           const int thresh_method, /// 0-Jurges, 1-Hopit
                           const int thresh_func, /// 0-exp, 1-lin
                           const int use_alpha, /// 0-auto
                           const double alpha_0,
                           const int negative, ///1 -yes 0 no
                           const int use_weights,
                           const Eigen::VectorXd weights){ /// 1-yes, 0-no

  Eigen::VectorXd reg_par = subvec(0, parcount(0) - 1L, parameters);
  int p01 = parcount(0) + parcount(1) - 1L;
  Eigen::VectorXd thresh_lambda = subvec(parcount(0), p01, parameters);
  Eigen::VectorXd thresh_gamma = subvec(p01 + 1L, p01 + parcount(2), parameters);
  Eigen::MatrixXd a = getThresholds(thresh_mm, thresh_lambda, thresh_gamma, thresh_no_cov,
                                    thresh_method, thresh_func, use_alpha, alpha_0 );
  Eigen::VectorXd b = reg_mm * reg_par;

  Eigen::VectorXd LO = colpath(a, yi, 0L) - b;
  Eigen::VectorXd HI = colpath(a, yi, 1L) - b;
  LO = exchange(LO, R_NegInf, -20);
  HI = exchange(HI, R_PosInf, 20);

  Eigen::VectorXd P;
  if (link == 0) {
    P = pstdnorm(HI) - pstdnorm(LO);
  } else {
    P = pstdlogit(HI) - pstdlogit(LO);
  }
  int d;
  if (negative == 1L) {
    d = -1;
  }  else d = 1;

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
                           const Eigen::MatrixXd reg_mm,
                           const Eigen::MatrixXd thresh_mm,
                           const Eigen::MatrixXd thresh_extd,
                           const Eigen::VectorXi parcount,
                           const int link, /// 0-probit , 1 - logit
                           const int thresh_no_cov, /// 0-NO , 1-Yes
                           const int thresh_method, /// 0-Jurges, 1-Hopit
                           const int thresh_func, /// 0-exp, 1-lin
                           const int use_alpha, /// 0-auto
                           const double alpha_0,
                           const int negative, ///1 -yes 0 no
                           const int use_weights, ///1 -yes 0 no
                           const Eigen::VectorXd weights){ /// 1-yes, 0-no

  long N = yi.size();

  Eigen::VectorXd reg_par = subvec(0, parcount(0) - 1L, parameters);
  int p01 = parcount(0) + parcount(1) - 1L;
  Eigen::VectorXd thresh_lambda = subvec(parcount(0), p01, parameters);
  Eigen::VectorXd thresh_gamma = subvec(p01 + 1L, p01 + parcount(2), parameters);
  Eigen::MatrixXd a = getThresholds(thresh_mm, thresh_lambda, thresh_gamma, thresh_no_cov,
                                    thresh_method, thresh_func, use_alpha, alpha_0 );
  Eigen::VectorXd b = reg_mm * reg_par;

  Eigen::VectorXd A2 = colpath(a, yi, 0L) - b;
  Eigen::VectorXd A1 = colpath(a, yi, 1L) - b;
  A2 = exchange(A2, R_NegInf, -20);
  A1 = exchange(A1, R_PosInf, 20);

  Eigen::VectorXd D, P, P1, P2, D1, D2;
  if (link == 0) {
    P1 = pstdnorm(A1);
    P2 = pstdnorm(A2);
    D1 = dstdnorm(A1);
    D2 = dstdnorm(A2);
  } else {
    P1 = pstdlogit(A1);
    P2 = pstdlogit(A2);
    D1 = dstdlogit(A1);
    D2 = dstdlogit(A2);
  }

  P = P1 - P2;
  D = D1 - D2;

  Eigen::VectorXd dlnLL_dX = P.array().cwiseInverse().matrix();

  Eigen::MatrixXd dlnLL_dbeta = Eigen::MatrixXd::Zero(D.size(), reg_par.size());
  Eigen::MatrixXd dd = - D.cwiseProduct(dlnLL_dX);
  for (int i = 0L; i < reg_par.size(); ++i) {
    dlnLL_dbeta.col(i) = dd;
  }

  dlnLL_dbeta = dlnLL_dbeta.cwiseProduct(reg_mm);

  Eigen::MatrixXd da = Eigen::MatrixXd::Ones(N, YYY1.cols());
  if (thresh_func == 0) {
    if (thresh_method==0) {
      da.col(0) = Eigen::VectorXd::Ones(N);
    } else {
      da.col(0) = a.col(1);
    }
    for (int i = 1; i < (a.cols() - 2); ++i){
      da.col(i) = a.col(i + 1) - a.col(i);
    }
  }

  Eigen::MatrixXd D1_(D1.rows(), da.cols()), D2_(D1.rows(), da.cols()), D3_(D1.rows(), da.cols());
  for (int i = 0; i < da.cols(); ++i){
    D1_.col(i) = D1;
    D2_.col(i) = D2;
    D3_.col(i) = dlnLL_dX;
  }
  Eigen::MatrixXd dlnLL_Lambda = (D1_.cwiseProduct(da.cwiseProduct(YYY1)) - D2_.cwiseProduct(da.cwiseProduct(YYY2))).cwiseProduct(D3_);

  Eigen::VectorXd res;
  if (use_weights==0){
    if (thresh_no_cov == 1L) {
      res.resize(dlnLL_dbeta.cols()+dlnLL_Lambda.cols());
      res << colSums_c(dlnLL_dbeta), colSums_c(dlnLL_Lambda);
    } else {
      Eigen::MatrixXd dlnLL_Gamma = rep_col_c(dlnLL_Lambda, thresh_mm.cols()).cwiseProduct(thresh_extd);
      res.resize(dlnLL_dbeta.cols()+dlnLL_Lambda.cols()+dlnLL_Gamma.cols());
      res << colSums_c(dlnLL_dbeta), colSums_c(dlnLL_Lambda), colSums_c(dlnLL_Gamma);
    }
  } else {
    if (thresh_no_cov == 1L) {
      res.resize(dlnLL_dbeta.cols()+dlnLL_Lambda.cols());
      res << colSums_weighted(dlnLL_dbeta, weights), colSums_weighted(dlnLL_Lambda, weights); // weigth sums by weights
    } else {
      Eigen::MatrixXd dlnLL_Gamma = rep_col_c(dlnLL_Lambda, thresh_mm.cols()).cwiseProduct(thresh_extd);
      res.resize(dlnLL_dbeta.cols()+dlnLL_Lambda.cols()+dlnLL_Gamma.cols());
      res << colSums_weighted(dlnLL_dbeta, weights), colSums_weighted(dlnLL_Lambda, weights), colSums_weighted(dlnLL_Gamma, weights);
    }
  }
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
                               const Eigen::MatrixXd reg_mm,
                               const Eigen::MatrixXd thresh_mm,
                               const Eigen::MatrixXd thresh_extd,
                               const Eigen::VectorXi parcount,
                               const int link, /// 0-probit , 1 - logit
                               const int thresh_no_cov, /// 0-NO , 1-Yes
                               const int thresh_method, /// 0-Jurges, 1-Hopit
                               const int thresh_func, /// 0-exp, 1-lin
                               const int use_alpha, /// 0-auto
                               const double alpha_0,
                               const int negative, ///1 -yes 0 no
                               const int use_weights, ///1 -yes 0 no
                               const Eigen::VectorXd weights){ /// 1-yes, 0-no

  long N = yi.size();

  Eigen::VectorXd reg_par = subvec(0, parcount(0) - 1L, parameters);
  int p01 = parcount(0) + parcount(1) - 1L;
  Eigen::VectorXd thresh_lambda = subvec(parcount(0), p01, parameters);
  Eigen::VectorXd thresh_gamma = subvec(p01 + 1L, p01 + parcount(2), parameters);
  Eigen::MatrixXd a = getThresholds(thresh_mm, thresh_lambda, thresh_gamma, thresh_no_cov,
                                    thresh_method, thresh_func, use_alpha, alpha_0 );
  Eigen::VectorXd b = reg_mm * reg_par;

  Eigen::VectorXd A2 = colpath(a, yi, 0L) - b;
  Eigen::VectorXd A1 = colpath(a, yi, 1L) - b;
  A2 = exchange(A2, R_NegInf, -20);
  A1 = exchange(A1, R_PosInf, 20);

  Eigen::VectorXd D, P, P1, P2, D1, D2;
  if (link == 0) {
    P1 = pstdnorm(A1);
    P2 = pstdnorm(A2);
    D1 = dstdnorm(A1);
    D2 = dstdnorm(A2);
  } else {
    P1 = pstdlogit(A1);
    P2 = pstdlogit(A2);
    D1 = dstdlogit(A1);
    D2 = dstdlogit(A2);
  }

  P = P1 - P2;
  D = D1 - D2;

  Eigen::VectorXd dlnLL_dX = P.array().cwiseInverse().matrix();

  Eigen::MatrixXd dlnLL_dbeta = Eigen::MatrixXd::Zero(D.size(), reg_par.size());
  Eigen::MatrixXd dd = - D.cwiseProduct(dlnLL_dX);
  for (int i = 0L; i < reg_par.size(); ++i) {
    dlnLL_dbeta.col(i) = dd;
  }

  dlnLL_dbeta = dlnLL_dbeta.cwiseProduct(reg_mm);

  Eigen::MatrixXd da = Eigen::MatrixXd::Ones(N, YYY1.cols());
  if (thresh_func == 0) {
    if (thresh_method==0) {
      da.col(0) = Eigen::VectorXd::Ones(N);
    } else {
      da.col(0) = a.col(1);
    }
    for (int i = 1; i < (a.cols() - 2); ++i){
      da.col(i) = a.col(i + 1) - a.col(i);
    }
  }

  Eigen::MatrixXd D1_(D1.rows(), da.cols()), D2_(D1.rows(), da.cols()), D3_(D1.rows(), da.cols());
  for (int i = 0; i < da.cols(); ++i){
    D1_.col(i) = D1;
    D2_.col(i) = D2;
    D3_.col(i) = dlnLL_dX;
  }

  Eigen::MatrixXd dlnLL_Lambda = (D1_.cwiseProduct(da.cwiseProduct(YYY1)) - D2_.cwiseProduct(da.cwiseProduct(YYY2))).cwiseProduct(D3_);

  int k;
  if (use_weights==0) k = 0; else k = parcount(2);
  Eigen::MatrixXd res(dlnLL_dbeta.rows(),dlnLL_dbeta.cols()+dlnLL_Lambda.cols()+k);

  if (use_weights==0){
    if (thresh_no_cov == 1L) {
      //res.resize(dlnLL_dbeta.rows(), dlnLL_dbeta.cols()+dlnLL_Lambda.cols());
      res << dlnLL_dbeta, dlnLL_Lambda;
    } else {
      Eigen::MatrixXd dlnLL_Gamma = rep_col_c(dlnLL_Lambda, thresh_mm.cols()).cwiseProduct(thresh_extd);
      //res.resize(dlnLL_dbeta.rows(),dlnLL_dbeta.cols()+dlnLL_Lambda.cols()+dlnLL_Gamma.cols());
      res << dlnLL_dbeta, dlnLL_Lambda, dlnLL_Gamma;
    }
  } else {
    if (thresh_no_cov == 1L) {
      //res.resize(dlnLL_dbeta.rows(),dlnLL_dbeta.cols()+dlnLL_Lambda.cols());
      res << weight_rows(dlnLL_dbeta, weights), weight_rows(dlnLL_Lambda, weights); // weigth sums by weights
    } else {
      Eigen::MatrixXd dlnLL_Gamma = rep_col_c(dlnLL_Lambda, thresh_mm.cols()).cwiseProduct(thresh_extd);
      //res.resize(dlnLL_dbeta.rows(),dlnLL_dbeta.cols()+dlnLL_Lambda.cols()+dlnLL_Gamma.cols());
      res << weight_rows(dlnLL_dbeta, weights), weight_rows(dlnLL_Lambda, weights), weight_rows(dlnLL_Gamma, weights);
    }
  }
  return(res);
}
