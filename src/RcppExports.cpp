// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// thread_affMtx
double thread_affMtx(unsigned int z_ini, unsigned int z_end, SEXP sexpX, bool isDistance, bool isSparse, SEXP sexpB, SEXP sexpP, SEXP sexpW);
RcppExport SEXP _bigMap_thread_affMtx(SEXP z_iniSEXP, SEXP z_endSEXP, SEXP sexpXSEXP, SEXP isDistanceSEXP, SEXP isSparseSEXP, SEXP sexpBSEXP, SEXP sexpPSEXP, SEXP sexpWSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type z_ini(z_iniSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type z_end(z_endSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sexpX(sexpXSEXP);
    Rcpp::traits::input_parameter< bool >::type isDistance(isDistanceSEXP);
    Rcpp::traits::input_parameter< bool >::type isSparse(isSparseSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sexpB(sexpBSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sexpP(sexpPSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sexpW(sexpWSEXP);
    rcpp_result_gen = Rcpp::wrap(thread_affMtx(z_ini, z_end, sexpX, isDistance, isSparse, sexpB, sexpP, sexpW));
    return rcpp_result_gen;
END_RCPP
}
// thread_repF
double thread_repF(unsigned int z_ini, unsigned int z_end, SEXP sexpY, double theta, SEXP sexpR);
RcppExport SEXP _bigMap_thread_repF(SEXP z_iniSEXP, SEXP z_endSEXP, SEXP sexpYSEXP, SEXP thetaSEXP, SEXP sexpRSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type z_ini(z_iniSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type z_end(z_endSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sexpY(sexpYSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sexpR(sexpRSEXP);
    rcpp_result_gen = Rcpp::wrap(thread_repF(z_ini, z_end, sexpY, theta, sexpR));
    return rcpp_result_gen;
END_RCPP
}
// thread_mIter
Rcpp::List thread_mIter(unsigned int z_ini, unsigned int z_end, SEXP sexpP, SEXP sexpW, SEXP sexpY, double sumP, double sumQ, SEXP sexpR, SEXP sexpU, SEXP sexpG, double eta, double alpha, double gain);
RcppExport SEXP _bigMap_thread_mIter(SEXP z_iniSEXP, SEXP z_endSEXP, SEXP sexpPSEXP, SEXP sexpWSEXP, SEXP sexpYSEXP, SEXP sumPSEXP, SEXP sumQSEXP, SEXP sexpRSEXP, SEXP sexpUSEXP, SEXP sexpGSEXP, SEXP etaSEXP, SEXP alphaSEXP, SEXP gainSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type z_ini(z_iniSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type z_end(z_endSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sexpP(sexpPSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sexpW(sexpWSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sexpY(sexpYSEXP);
    Rcpp::traits::input_parameter< double >::type sumP(sumPSEXP);
    Rcpp::traits::input_parameter< double >::type sumQ(sumQSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sexpR(sexpRSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sexpU(sexpUSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sexpG(sexpGSEXP);
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type gain(gainSEXP);
    rcpp_result_gen = Rcpp::wrap(thread_mIter(z_ini, z_end, sexpP, sexpW, sexpY, sumP, sumQ, sexpR, sexpU, sexpG, eta, alpha, gain));
    return rcpp_result_gen;
END_RCPP
}
// z_kNP
arma::Mat<int> z_kNP(int thread_rank, int threads, SEXP sexpX, SEXP sexpY, bool is_distance, bool is_sparse, const arma::Col<int>& K, double sampling);
RcppExport SEXP _bigMap_z_kNP(SEXP thread_rankSEXP, SEXP threadsSEXP, SEXP sexpXSEXP, SEXP sexpYSEXP, SEXP is_distanceSEXP, SEXP is_sparseSEXP, SEXP KSEXP, SEXP samplingSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type thread_rank(thread_rankSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sexpX(sexpXSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sexpY(sexpYSEXP);
    Rcpp::traits::input_parameter< bool >::type is_distance(is_distanceSEXP);
    Rcpp::traits::input_parameter< bool >::type is_sparse(is_sparseSEXP);
    Rcpp::traits::input_parameter< const arma::Col<int>& >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type sampling(samplingSEXP);
    rcpp_result_gen = Rcpp::wrap(z_kNP(thread_rank, threads, sexpX, sexpY, is_distance, is_sparse, K, sampling));
    return rcpp_result_gen;
END_RCPP
}
// z_hlCorr
double z_hlCorr(SEXP sexpX, SEXP sexpY, int zSampleSize, bool is_distance, bool is_sparse);
RcppExport SEXP _bigMap_z_hlCorr(SEXP sexpXSEXP, SEXP sexpYSEXP, SEXP zSampleSizeSEXP, SEXP is_distanceSEXP, SEXP is_sparseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type sexpX(sexpXSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sexpY(sexpYSEXP);
    Rcpp::traits::input_parameter< int >::type zSampleSize(zSampleSizeSEXP);
    Rcpp::traits::input_parameter< bool >::type is_distance(is_distanceSEXP);
    Rcpp::traits::input_parameter< bool >::type is_sparse(is_sparseSEXP);
    rcpp_result_gen = Rcpp::wrap(z_hlCorr(sexpX, sexpY, zSampleSize, is_distance, is_sparse));
    return rcpp_result_gen;
END_RCPP
}
// grid_init
arma::Mat<double> grid_init(arma::Col<double> X, arma::Col<double> Y);
RcppExport SEXP _bigMap_grid_init(SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Col<double> >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::Col<double> >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(grid_init(X, Y));
    return rcpp_result_gen;
END_RCPP
}
// grid_p2cell
arma::Col<int> grid_p2cell(double x, double y, arma::Mat<double> grid);
RcppExport SEXP _bigMap_grid_p2cell(SEXP xSEXP, SEXP ySEXP, SEXP gridSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::Mat<double> >::type grid(gridSEXP);
    rcpp_result_gen = Rcpp::wrap(grid_p2cell(x, y, grid));
    return rcpp_result_gen;
END_RCPP
}
// grid_D2cell
arma::Mat<double> grid_D2cell(arma::Mat<double> D, arma::Mat<double> grid);
RcppExport SEXP _bigMap_grid_D2cell(SEXP DSEXP, SEXP gridSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Mat<double> >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::Mat<double> >::type grid(gridSEXP);
    rcpp_result_gen = Rcpp::wrap(grid_D2cell(D, grid));
    return rcpp_result_gen;
END_RCPP
}
// grid_n2cell
arma::Col<int> grid_n2cell(int n, arma::Mat<double> grid);
RcppExport SEXP _bigMap_grid_n2cell(SEXP nSEXP, SEXP gridSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::Mat<double> >::type grid(gridSEXP);
    rcpp_result_gen = Rcpp::wrap(grid_n2cell(n, grid));
    return rcpp_result_gen;
END_RCPP
}
// grid_N2cell
arma::Mat<int> grid_N2cell(arma::Mat<double> grid);
RcppExport SEXP _bigMap_grid_N2cell(SEXP gridSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Mat<double> >::type grid(gridSEXP);
    rcpp_result_gen = Rcpp::wrap(grid_N2cell(grid));
    return rcpp_result_gen;
END_RCPP
}
// grid_M2cell
arma::Mat<int> grid_M2cell(arma::Col<int> M, arma::Mat<double> grid);
RcppExport SEXP _bigMap_grid_M2cell(SEXP MSEXP, SEXP gridSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Col<int> >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::Mat<double> >::type grid(gridSEXP);
    rcpp_result_gen = Rcpp::wrap(grid_M2cell(M, grid));
    return rcpp_result_gen;
END_RCPP
}
// grid_bound
arma::Col<int> grid_bound(int n, arma::Mat<double> grid);
RcppExport SEXP _bigMap_grid_bound(SEXP nSEXP, SEXP gridSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::Mat<double> >::type grid(gridSEXP);
    rcpp_result_gen = Rcpp::wrap(grid_bound(n, grid));
    return rcpp_result_gen;
END_RCPP
}
// grid_cross
arma::Col<int> grid_cross(int n, arma::Mat<double> grid);
RcppExport SEXP _bigMap_grid_cross(SEXP nSEXP, SEXP gridSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::Mat<double> >::type grid(gridSEXP);
    rcpp_result_gen = Rcpp::wrap(grid_cross(n, grid));
    return rcpp_result_gen;
END_RCPP
}
// grid_peaks
arma::Col<int> grid_peaks(arma::Mat<double> Z, arma::Mat<double> grid);
RcppExport SEXP _bigMap_grid_peaks(SEXP ZSEXP, SEXP gridSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Mat<double> >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< arma::Mat<double> >::type grid(gridSEXP);
    rcpp_result_gen = Rcpp::wrap(grid_peaks(Z, grid));
    return rcpp_result_gen;
END_RCPP
}
// wtt_cpp
Rcpp::List wtt_cpp(arma::Col<double> X, arma::Col<double> Y, arma::Mat<double> Z);
RcppExport SEXP _bigMap_wtt_cpp(SEXP XSEXP, SEXP YSEXP, SEXP ZSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Col<double> >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::Col<double> >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::Mat<double> >::type Z(ZSEXP);
    rcpp_result_gen = Rcpp::wrap(wtt_cpp(X, Y, Z));
    return rcpp_result_gen;
END_RCPP
}
// zBeta
Rcpp::NumericMatrix zBeta(int thread_rank, int threads, SEXP sexpX, bool is_distance, bool is_sparse, int ppx, double xppx);
RcppExport SEXP _bigMap_zBeta(SEXP thread_rankSEXP, SEXP threadsSEXP, SEXP sexpXSEXP, SEXP is_distanceSEXP, SEXP is_sparseSEXP, SEXP ppxSEXP, SEXP xppxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type thread_rank(thread_rankSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sexpX(sexpXSEXP);
    Rcpp::traits::input_parameter< bool >::type is_distance(is_distanceSEXP);
    Rcpp::traits::input_parameter< bool >::type is_sparse(is_sparseSEXP);
    Rcpp::traits::input_parameter< int >::type ppx(ppxSEXP);
    Rcpp::traits::input_parameter< double >::type xppx(xppxSEXP);
    rcpp_result_gen = Rcpp::wrap(zBeta(thread_rank, threads, sexpX, is_distance, is_sparse, ppx, xppx));
    return rcpp_result_gen;
END_RCPP
}
// centerScale
void centerScale(SEXP sexpX, bool is_distance, bool is_sparse);
RcppExport SEXP _bigMap_centerScale(SEXP sexpXSEXP, SEXP is_distanceSEXP, SEXP is_sparseSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type sexpX(sexpXSEXP);
    Rcpp::traits::input_parameter< bool >::type is_distance(is_distanceSEXP);
    Rcpp::traits::input_parameter< bool >::type is_sparse(is_sparseSEXP);
    centerScale(sexpX, is_distance, is_sparse);
    return R_NilValue;
END_RCPP
}
// sckt_zTSNE
double sckt_zTSNE(unsigned int thread_rank, unsigned int threads, unsigned int layers, SEXP sexpX, SEXP sexpB, SEXP sexpY, SEXP sexpI, int iters, double nnSize, double theta, double lRate, double alpha, double gain, bool isDistance, bool isSparse);
RcppExport SEXP _bigMap_sckt_zTSNE(SEXP thread_rankSEXP, SEXP threadsSEXP, SEXP layersSEXP, SEXP sexpXSEXP, SEXP sexpBSEXP, SEXP sexpYSEXP, SEXP sexpISEXP, SEXP itersSEXP, SEXP nnSizeSEXP, SEXP thetaSEXP, SEXP lRateSEXP, SEXP alphaSEXP, SEXP gainSEXP, SEXP isDistanceSEXP, SEXP isSparseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type thread_rank(thread_rankSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type layers(layersSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sexpX(sexpXSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sexpB(sexpBSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sexpY(sexpYSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sexpI(sexpISEXP);
    Rcpp::traits::input_parameter< int >::type iters(itersSEXP);
    Rcpp::traits::input_parameter< double >::type nnSize(nnSizeSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type lRate(lRateSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type gain(gainSEXP);
    Rcpp::traits::input_parameter< bool >::type isDistance(isDistanceSEXP);
    Rcpp::traits::input_parameter< bool >::type isSparse(isSparseSEXP);
    rcpp_result_gen = Rcpp::wrap(sckt_zTSNE(thread_rank, threads, layers, sexpX, sexpB, sexpY, sexpI, iters, nnSize, theta, lRate, alpha, gain, isDistance, isSparse));
    return rcpp_result_gen;
END_RCPP
}
// zChnks
void zChnks(Rcpp::List& Z_list, const arma::Mat<double>& Y, const arma::Col<int>& I, const Rcpp::List& brks_list);
RcppExport SEXP _bigMap_zChnks(SEXP Z_listSEXP, SEXP YSEXP, SEXP ISEXP, SEXP brks_listSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List& >::type Z_list(Z_listSEXP);
    Rcpp::traits::input_parameter< const arma::Mat<double>& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::Col<int>& >::type I(ISEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type brks_list(brks_listSEXP);
    zChnks(Z_list, Y, I, brks_list);
    return R_NilValue;
END_RCPP
}
// updateY
void updateY(arma::Mat<double>& Y, const arma::Col<int>& I, const Rcpp::List& zMap_list, const Rcpp::List& brks_list);
RcppExport SEXP _bigMap_updateY(SEXP YSEXP, SEXP ISEXP, SEXP zMap_listSEXP, SEXP brks_listSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Mat<double>& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::Col<int>& >::type I(ISEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type zMap_list(zMap_listSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type brks_list(brks_listSEXP);
    updateY(Y, I, zMap_list, brks_list);
    return R_NilValue;
END_RCPP
}
// mpi_zTSNE
double mpi_zTSNE(unsigned int thread_rank, SEXP sexpX, SEXP sexpB, arma::Mat<double>& Y, arma::Col<int> indexes, int iters, double nnSize, double theta, double lRate, double alpha, double gain, bool isDistance, bool isSparse);
RcppExport SEXP _bigMap_mpi_zTSNE(SEXP thread_rankSEXP, SEXP sexpXSEXP, SEXP sexpBSEXP, SEXP YSEXP, SEXP indexesSEXP, SEXP itersSEXP, SEXP nnSizeSEXP, SEXP thetaSEXP, SEXP lRateSEXP, SEXP alphaSEXP, SEXP gainSEXP, SEXP isDistanceSEXP, SEXP isSparseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type thread_rank(thread_rankSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sexpX(sexpXSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sexpB(sexpBSEXP);
    Rcpp::traits::input_parameter< arma::Mat<double>& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::Col<int> >::type indexes(indexesSEXP);
    Rcpp::traits::input_parameter< int >::type iters(itersSEXP);
    Rcpp::traits::input_parameter< double >::type nnSize(nnSizeSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type lRate(lRateSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type gain(gainSEXP);
    Rcpp::traits::input_parameter< bool >::type isDistance(isDistanceSEXP);
    Rcpp::traits::input_parameter< bool >::type isSparse(isSparseSEXP);
    rcpp_result_gen = Rcpp::wrap(mpi_zTSNE(thread_rank, sexpX, sexpB, Y, indexes, iters, nnSize, theta, lRate, alpha, gain, isDistance, isSparse));
    return rcpp_result_gen;
END_RCPP
}
// nnSS_chk
Rcpp::List nnSS_chk(SEXP sexpX, SEXP sexpB, arma::Col<int> indexes, bool isDistance, bool isSparse, unsigned int nnSize);
RcppExport SEXP _bigMap_nnSS_chk(SEXP sexpXSEXP, SEXP sexpBSEXP, SEXP indexesSEXP, SEXP isDistanceSEXP, SEXP isSparseSEXP, SEXP nnSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type sexpX(sexpXSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sexpB(sexpBSEXP);
    Rcpp::traits::input_parameter< arma::Col<int> >::type indexes(indexesSEXP);
    Rcpp::traits::input_parameter< bool >::type isDistance(isDistanceSEXP);
    Rcpp::traits::input_parameter< bool >::type isSparse(isSparseSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nnSize(nnSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(nnSS_chk(sexpX, sexpB, indexes, isDistance, isSparse, nnSize));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bigMap_thread_affMtx", (DL_FUNC) &_bigMap_thread_affMtx, 8},
    {"_bigMap_thread_repF", (DL_FUNC) &_bigMap_thread_repF, 5},
    {"_bigMap_thread_mIter", (DL_FUNC) &_bigMap_thread_mIter, 13},
    {"_bigMap_z_kNP", (DL_FUNC) &_bigMap_z_kNP, 8},
    {"_bigMap_z_hlCorr", (DL_FUNC) &_bigMap_z_hlCorr, 5},
    {"_bigMap_grid_init", (DL_FUNC) &_bigMap_grid_init, 2},
    {"_bigMap_grid_p2cell", (DL_FUNC) &_bigMap_grid_p2cell, 3},
    {"_bigMap_grid_D2cell", (DL_FUNC) &_bigMap_grid_D2cell, 2},
    {"_bigMap_grid_n2cell", (DL_FUNC) &_bigMap_grid_n2cell, 2},
    {"_bigMap_grid_N2cell", (DL_FUNC) &_bigMap_grid_N2cell, 1},
    {"_bigMap_grid_M2cell", (DL_FUNC) &_bigMap_grid_M2cell, 2},
    {"_bigMap_grid_bound", (DL_FUNC) &_bigMap_grid_bound, 2},
    {"_bigMap_grid_cross", (DL_FUNC) &_bigMap_grid_cross, 2},
    {"_bigMap_grid_peaks", (DL_FUNC) &_bigMap_grid_peaks, 2},
    {"_bigMap_wtt_cpp", (DL_FUNC) &_bigMap_wtt_cpp, 3},
    {"_bigMap_zBeta", (DL_FUNC) &_bigMap_zBeta, 7},
    {"_bigMap_centerScale", (DL_FUNC) &_bigMap_centerScale, 3},
    {"_bigMap_sckt_zTSNE", (DL_FUNC) &_bigMap_sckt_zTSNE, 15},
    {"_bigMap_zChnks", (DL_FUNC) &_bigMap_zChnks, 4},
    {"_bigMap_updateY", (DL_FUNC) &_bigMap_updateY, 4},
    {"_bigMap_mpi_zTSNE", (DL_FUNC) &_bigMap_mpi_zTSNE, 13},
    {"_bigMap_nnSS_chk", (DL_FUNC) &_bigMap_nnSS_chk, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_bigMap(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
