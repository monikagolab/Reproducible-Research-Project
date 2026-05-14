"""Data preparation and stationarity testing for cryptocurrency returns."""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import yfinance as yf
from arch.unitroot import PhillipsPerron
from statsmodels.tsa.stattools import adfuller, kpss


class CryptoDataLoader:
    """Download cryptocurrency price data from Yahoo Finance."""

    def __init__(
        self,
        tickers: list[str] | None = None,
        start_date: str = "2017-01-01",
        end_date: str = "2025-09-06",
    ) -> None:
        """Initialize the data loader.

        Args:
            tickers: Yahoo Finance tickers.
            start_date: Start date for the sample.
            end_date: End date for the sample.
        """
        self.tickers = tickers or ["BTC-USD", "ETH-USD", "SOL-USD"]
        self.start_date = start_date
        self.end_date = end_date

    def load_prices(self) -> pd.DataFrame:
        """Download daily closing prices and rename columns.

        Returns:
            DataFrame with BTC, ETH and SOL prices.
        """
        prices = yf.download(
            self.tickers,
            start=self.start_date,
            end=self.end_date,
            auto_adjust=False,
        )["Close"]

        prices.columns = ["BTC", "ETH", "SOL"]
        return prices


class CryptoPreprocessor:
    """Prepare cryptocurrency price data for time-series modelling."""

    @staticmethod
    def compute_log_returns(prices: pd.DataFrame) -> pd.DataFrame:
        """Compute daily log returns and remove missing values.

        Args:
            prices: DataFrame with cryptocurrency prices.

        Returns:
            DataFrame with aligned daily log returns.
        """
        returns = np.log(prices / prices.shift(1)).dropna()
        return returns


class StationarityAnalyzer:
    """Run ADF, PP and KPSS stationarity tests."""

    @staticmethod
    def adf_test(series: pd.Series) -> tuple[float, float]:
        """Run the Augmented Dickey-Fuller test.

        Args:
            series: Time series to test.

        Returns:
            Test statistic and p-value.
        """
        result = adfuller(series, regression="c", autolag="AIC")
        return result[0], result[1]

    @staticmethod
    def pp_test(series: pd.Series) -> tuple[float, float]:
        """Run the Phillips-Perron test.

        Args:
            series: Time series to test.

        Returns:
            Test statistic and p-value.
        """
        result = PhillipsPerron(series)
        return result.stat, result.pvalue

    @staticmethod
    def kpss_test(series: pd.Series) -> tuple[float, float]:
        """Run the KPSS test.

        Args:
            series: Time series to test.

        Returns:
            Test statistic and p-value.
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = kpss(series, regression="c", nlags="auto")

        return result[0], result[1]

    def run_all_tests(self, returns: pd.DataFrame) -> pd.DataFrame:
        """Run ADF, PP and KPSS tests for BTC, ETH and SOL.

        Args:
            returns: DataFrame with daily log returns.

        Returns:
            DataFrame with test statistics and p-values.
        """
        results = []

        for col in ["BTC", "ETH", "SOL"]:
            adf_stat, adf_p = self.adf_test(returns[col])
            pp_stat, pp_p = self.pp_test(returns[col])
            kpss_stat, kpss_p = self.kpss_test(returns[col])

            results.append([col, "ADF", adf_stat, adf_p])
            results.append([col, "PP", pp_stat, pp_p])
            results.append([col, "KPSS", kpss_stat, kpss_p])

        return pd.DataFrame(
            results,
            columns=["Variable", "Test", "Statistic", "p-value"],
        )

    @staticmethod
    def summarize_results(results_df: pd.DataFrame) -> pd.DataFrame:
        """Create a compact stationarity summary table.

        Args:
            results_df: Full stationarity test results.

        Returns:
            Summary table with final stationarity conclusion.
        """
        summary_rows = []

        for col in ["BTC", "ETH", "SOL"]:
            sub = results_df[results_df["Variable"] == col]

            adf_p = sub[sub["Test"] == "ADF"]["p-value"].values[0]
            pp_p = sub[sub["Test"] == "PP"]["p-value"].values[0]
            kpss_p = sub[sub["Test"] == "KPSS"]["p-value"].values[0]

            conclusion = (
                "Stationary"
                if adf_p < 0.05 and pp_p < 0.05 and kpss_p > 0.05
                else "Needs review"
            )

            summary_rows.append([col, adf_p, pp_p, kpss_p, conclusion])

        return pd.DataFrame(
            summary_rows,
            columns=[
                "Variable",
                "ADF p-value",
                "PP p-value",
                "KPSS p-value",
                "Conclusion",
            ],
        )


class DataStationarityPipeline:
    """Full pipeline for data loading, preprocessing and stationarity testing."""

    def __init__(
        self,
        tickers: list[str] | None = None,
        start_date: str = "2017-01-01",
        end_date: str = "2025-09-06",
    ) -> None:
        """Initialize the full data and stationarity pipeline."""
        self.loader = CryptoDataLoader(tickers, start_date, end_date)
        self.preprocessor = CryptoPreprocessor()
        self.analyzer = StationarityAnalyzer()

    def run(self) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """Run the complete pipeline.

        Returns:
            Prices, returns, full test results and summary table.
        """
        prices = self.loader.load_prices()
        returns = self.preprocessor.compute_log_returns(prices)
        results = self.analyzer.run_all_tests(returns)
        summary = self.analyzer.summarize_results(results)

        return prices, returns, results, summary