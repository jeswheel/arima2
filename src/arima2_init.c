#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP ARIMA_CSS(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ARIMA_Gradtrans(SEXP, SEXP);
extern SEXP ARIMA_Invtrans(SEXP, SEXP);
extern SEXP ARIMA_Like(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ARIMA_transPars(SEXP, SEXP, SEXP);
extern SEXP ARIMA_undoPars(SEXP, SEXP);
extern SEXP getQ0(SEXP, SEXP);
extern SEXP getQ0bis(SEXP, SEXP, SEXP);
extern SEXP TSconv(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"ARIMA_CSS",       (DL_FUNC) &ARIMA_CSS,       6},
    {"ARIMA_Gradtrans", (DL_FUNC) &ARIMA_Gradtrans, 2},
    {"ARIMA_Invtrans",  (DL_FUNC) &ARIMA_Invtrans,  2},
    {"ARIMA_Like",      (DL_FUNC) &ARIMA_Like,      5},
    {"ARIMA_transPars", (DL_FUNC) &ARIMA_transPars, 3},
    {"ARIMA_undoPars",  (DL_FUNC) &ARIMA_undoPars,  2},
    {"getQ0",           (DL_FUNC) &getQ0,           2},
    {"getQ0bis",        (DL_FUNC) &getQ0bis,        3},
    {"TSconv",          (DL_FUNC) &TSconv,          2},
    {NULL, NULL, 0}
};

void R_init_arima2(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
