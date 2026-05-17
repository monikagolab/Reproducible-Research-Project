"""Granger causality testing and out-of-sample forecast validation."""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
from statsmodels.tsa.api import VAR
from statsmodels.tsa.stattools import grangercausalitytests

logger = logging.getLogger(__name__)


class GrangerTester:
    """Tests pairwise Granger causality for all asset combinations.

    Replicates R: causality(var_model, cause=X)$Granger for each asset.

    Attributes:
        data: DataFrame of log-returns (T x k).
        max_lags: Maximum lag order for the Granger F-test.
    """

    def __init__(self, data: pd.DataFrame, max_lags: int = 5) -> None:
        """Initialize with log-return data.

        Args:
            data: DataFrame of log-returns, one column per asset.
            max_lags: Maximum lag order to test.
        """
        self.data = data.copy()
        self.max_lags = max_lags

    def test_pair(self, cause: str, effect: str) -> dict:
        """Test whether 'cause' Granger-causes 'effect'.

        Uses the F-test (ssr_ftest) at the lag with the lowest p-value.

        Args:
            cause: Column name of the potential causing variable.
            effect: Column name of the potentially affected variable.

        Returns:
            Dict with cause, effect, best_lag, f_stat, p_value, conclusion.
        """
        xy = self.data[[effect, cause]].dropna()
        results = grangercausalitytests(xy, maxlag=self.max_lags, verbose=False)
        # Pick the lag with the best (lowest) p-value
        best_lag, best_stats = min(
            results.items(),
            key=lambda item: item[1][0]["ssr_ftest"][1],
        )
        f_stat, p_val = best_stats[0]["ssr_ftest"][:2]
        return {
            "cause": cause,
            "effect": effect,
            "best_lag": best_lag,
            "f_stat": round(float(f_stat), 4),
            "p_value": round(float(p_val), 4),
            "conclusion": ("Granger-causes" if p_val < 0.05 else "No causality"),
        }

    def test_all(self) -> pd.DataFrame:
        """Run all ordered pairwise Granger tests.

        Returns:
            DataFrame with one row per (cause, effect) ordered pair.
        """
        cols = list(self.data.columns)
        rows = [
            self.test_pair(cause, effect)
            for cause in cols
            for effect in cols
            if cause != effect
        ]
        return pd.DataFrame(rows)


class OOSValidator:
    """Rolling expanding-window out-of-sample forecast evaluation.

    Compares three models:
        VAR(p)  - the main model
        AR(1)   - univariate autoregressive benchmark
        Zero    - naive zero-return forecast (random walk benchmark)

    Replicates R Section 6: split at 80%, then for each step t in
    [split, n-1]: fit on [0..t], forecast t+1, compare to actual.

    Attributes:
        data: Full log-return DataFrame (T x k).
        p: VAR lag order (should match the fitted VARAnalysis.p).
        split: Index of the last training observation.
    """

    def __init__(
        self,
        data: pd.DataFrame,
        p: int,
        train_fraction: float = 0.8,
    ) -> None:
        """Initialise with data, lag order, and train/test split.

        Args:
            data: Full log-return DataFrame.
            p: VAR lag order.
            train_fraction: Fraction of data in the initial training window.
        """
        self.data = data.copy()
        self.p = p
        self.split = max(p + 50, int(train_fraction * len(data)))
        self.columns = list(data.columns)

    @staticmethod
    def _mae(errors: np.ndarray) -> float:
        """Mean Absolute Error."""
        return float(np.mean(np.abs(errors)))

    @staticmethod
    def _rmse(errors: np.ndarray) -> float:
        """Root Mean Squared Error."""
        return float(np.sqrt(np.mean(errors**2)))

    def _var_forecast_1(self, train: pd.DataFrame) -> dict[str, float]:
        """One-step-ahead VAR forecast.

        Args:
            train: Training slice of data ending at current time t.

        Returns:
            Dict mapping column name to forecast value.
        """
        fit = VAR(train).fit(maxlags=self.p, ic=None, trend="c")
        fc = fit.forecast(train.values[-self.p :], steps=1)[0]
        return dict(zip(self.columns, fc, strict=False))

    def _ar1_forecast_1(self, y: np.ndarray) -> float:
        """One-step-ahead AR(1) forecast via OLS.

        Args:
            y: Full training history for one series up to time t.

        Returns:
            Scalar one-step-ahead forecast.
        """
        if len(y) < 5:
            return float(np.nan)
        b = np.polyfit(y[:-1], y[1:], deg=1)  # [slope, intercept]
        return float(b[0] * y[-1] + b[1])

    def run(self) -> pd.DataFrame:
        """Execute the full rolling OOS evaluation loop.

        Returns:
            DataFrame with columns: date, actual values, VAR/AR1/ZERO
            forecasts for each asset, and error metrics.
        """
        n = len(self.data)

        records = []
        for t_idx, t in enumerate(range(self.split, n - 1)):
            train = self.data.iloc[: t + 1]
            actual = self.data.iloc[t + 1]

            var_fc = self._var_forecast_1(train)
            row: dict = {"date": self.data.index[t + 1]}

            for col in self.columns:
                row[f"actual_{col}"] = float(actual[col])
                row[f"VAR_{col}"] = var_fc[col]
                row[f"AR1_{col}"] = self._ar1_forecast_1(train[col].values)
                row[f"ZERO_{col}"] = 0.0

            records.append(row)
            if (t_idx + 1) % 500 == 0:
                logger.info("OOS step %d / %d", t_idx + 1, n - self.split - 1)

        return pd.DataFrame(records)

    def metrics(self, oos_df: pd.DataFrame) -> pd.DataFrame:
        """Compute MAE and RMSE for each model and asset.

        Args:
            oos_df: DataFrame returned by run().

        Returns:
            DataFrame with columns: model, series, mae, rmse.
            Also includes a row for the mean across series per model.
        """
        models = ["VAR", "AR1", "ZERO"]
        rows = []
        for model in models:
            for col in self.columns:
                errors = (oos_df[f"actual_{col}"] - oos_df[f"{model}_{col}"]).values
                rows.append(
                    {
                        "model": model,
                        "series": col,
                        "mae": round(self._mae(errors), 6),
                        "rmse": round(self._rmse(errors), 6),
                    }
                )
        detail = pd.DataFrame(rows)
        avg = detail.groupby("model")[["mae", "rmse"]].mean().round(6).reset_index()
        avg["series"] = "AVERAGE"
        return pd.concat([detail, avg], ignore_index=True)
