"""Stationarity tests and VAR model estimation with diagnostics."""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
from statsmodels.stats.diagnostic import acorr_ljungbox
from statsmodels.tsa.api import VAR
from statsmodels.tsa.stattools import adfuller, kpss

logger = logging.getLogger(__name__)


class StationarityTester:
    """Runs ADF and KPSS unit-root tests on a univariate time series.

    Replicates R:
        ur.df(x, type='drift', selectlags='AIC')   -> ADF
        ur.kpss(x, type='mu', lags='short')         -> KPSS

    ADF null hypothesis  : unit root (non-stationary). Reject -> stationary.
    KPSS null hypothesis : stationary. Reject -> non-stationary.

    Attributes:
        series: The univariate series under test (NaNs dropped).
        name: Label used in printed output.
    """

    def __init__(self, series: pd.Series, name: str = "series") -> None:
        """Initialise with a pandas Series.

        Args:
            series: Univariate time series to test.
            name: Human-readable label for reports.
        """
        self.series = series.dropna()
        self.name = name

    def adf_test(self) -> dict:
        """Augmented Dickey-Fuller test with automatic AIC lag selection.

        Returns:
            Dict with test, series, statistic, p_value,
            lags_used, critical_values, conclusion.
        """
        stat, p_val, lags, _, crit, _ = adfuller(
            self.series, regression="c", autolag="AIC"
        )
        return {
            "test": "ADF",
            "series": self.name,
            "statistic": round(stat, 4),
            "p_value": round(p_val, 4),
            "lags_used": lags,
            "critical_values": {k: round(v, 4) for k, v in crit.items()},
            "conclusion": "Stationary" if p_val < 0.05 else "Non-stationary",
        }

    def kpss_test(self) -> dict:
        """KPSS test for level stationarity.

        Returns:
            Dict with test, series, statistic, p_value,
            lags_used, critical_values, conclusion.
        """
        stat, p_val, lags, crit = kpss(
            self.series, regression="c", nlags="auto"
        )
        return {
            "test": "KPSS",
            "series": self.name,
            "statistic": round(stat, 4),
            "p_value": round(p_val, 4),
            "lags_used": lags,
            "critical_values": {k: round(v, 4) for k, v in crit.items()},
            "conclusion": "Non-stationary" if p_val < 0.05 else "Stationary",
        }

    def summary(self) -> pd.DataFrame:
        """Run both tests and return a combined summary table.

        Returns:
            DataFrame with one row per test.
        """
        rows = [self.adf_test(), self.kpss_test()]
        return pd.DataFrame(rows)[
            ["test", "series", "statistic", "p_value", "conclusion"]
        ]


class VARAnalysis:
    """Fits a VAR(p) model with lag selection, stability check, and diagnostics.

    Replicates R workflow:
        VARselect -> VAR -> roots -> serial.test -> cor(residuals)
        -> irf -> fevd

    Lag order selected by BIC (SC in R notation), falling back to AIC
    if BIC returns a value less than 1.

    Attributes:
        data: DataFrame of log-returns (T x k).
        max_lags: Upper bound for lag-order search.
        p: Selected lag order (set after fit()).
        model: Fitted statsmodels VARResultsWrapper (set after fit()).
    """

    def __init__(self, data: pd.DataFrame, max_lags: int = 10) -> None:
        """Initialise with multivariate log-return data.

        Args:
            data: DataFrame where each column is one asset's log-returns.
            max_lags: Maximum lag order to evaluate during selection.
        """
        self.data = data.copy()
        self.max_lags = max_lags
        self.p: int | None = None
        self.model = None

    def fit(self) -> None:
        """Select lag order by BIC->AIC fallback and fit VAR(p) with constant.

        Sets self.p and self.model.
        Replicates R: VAR(ret3, p=p, type='const').
        """
        var = VAR(self.data)
        sel = var.select_order(maxlags=self.max_lags)
        logger.info("Lag selection results:\n%s", sel.summary())

        p_bic = int(sel.bic)
        p_aic = int(sel.aic)
        self.p = p_bic if p_bic >= 1 else max(p_aic, 1)

        self.model = var.fit(maxlags=self.p, ic=None, trend="c")
        logger.info(
            "VAR(%d) fitted | total obs: %d | used: %d",
            self.p,
            len(self.data),
            len(self.data) - self.p,
        )

    def _require_fit(self) -> None:
        """Raise an error if fit() has not been called yet."""
        if self.model is None:
            raise RuntimeError("Call fit() before using this method.")

    def is_stable(self) -> bool:
        """Check VAR stability: all companion matrix eigenvalue moduli < 1.

        Replicates R: max(abs(roots(var_model))) < 1.

        Returns:
            True if the model is stable, False otherwise.
        """
        self._require_fit()
        max_root = float(np.max(np.abs(self.model.roots)))
        stable = max_root < 1.0
        logger.info(
            "Max |root| = %.5f  ->  %s",
            max_root,
            "STABLE" if stable else "UNSTABLE",
        )
        return stable

    def residual_correlation(self) -> pd.DataFrame:
        """Correlation matrix of VAR residuals.

        Replicates R: cor(residuals(var_model)).

        Returns:
            k x k DataFrame of pairwise residual correlations.
        """
        self._require_fit()
        resid = pd.DataFrame(self.model.resid, columns=self.data.columns)
        return resid.corr().round(4)

    def portmanteau_test(self, lags: int = 10) -> pd.DataFrame:
        """Per-equation Ljung-Box test for residual autocorrelation.

        Approximates R: serial.test(var_model, lags.pt=10).

        Args:
            lags: Number of lags to include in the test statistic.

        Returns:
            DataFrame with columns: equation, lb_stat, lb_pvalue.
        """
        self._require_fit()
        rows = []
        for col in self.data.columns:
            lb = acorr_ljungbox(
                self.model.resid[col], lags=[lags], return_df=True
            )
            rows.append({
                "equation": col,
                "lb_stat": round(float(lb["lb_stat"].iloc[-1]), 4),
                "lb_pvalue": round(float(lb["lb_pvalue"].iloc[-1]), 4),
            })
        return pd.DataFrame(rows)

    def irf(self, periods: int = 3):
        """Compute Impulse Response Functions (Cholesky orthogonalised).

        Replicates R: irf(var_model, n.ahead=3, ortho=TRUE).

        Args:
            periods: Number of periods ahead for impulse responses.

        Returns:
            statsmodels IRAnalysis object. Call .plot() on the result.
        """
        self._require_fit()
        return self.model.irf(periods=periods)

    def fevd(self, periods: int = 3):
        """Compute Forecast Error Variance Decomposition.

        Replicates R: fevd(var_model, n.ahead=3).

        Args:
            periods: Forecast horizon for the decomposition.

        Returns:
            statsmodels FEVD object. Call .plot() on the result.
        """
        self._require_fit()
        return self.model.fevd(periods=periods)