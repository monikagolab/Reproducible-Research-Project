"""Cryptocurrency price data downloading and log-return computation."""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
import yfinance as yf

logger = logging.getLogger(__name__)


class DataLoader:
    """Downloads adjusted close prices and computes log-returns.

    Replicates the R script's data preparation:
    px <- Reduce(merge, lapply(tickers, get_px))
    ret3 <- na.omit(diff(log(px)))

    Attributes:
        tickers: Yahoo Finance ticker symbols.
        start_date: First date of data range (YYYY-MM-DD).
        end_date: Last date of data range (YYYY-MM-DD).
    """

    TICKER_MAP: dict[str, str] = {
        "BTC-USD": "BTC",
        "ETH-USD": "ETH",
        "SOL-USD": "SOL",
    }

    def __init__(
        self,
        tickers: list[str] | None = None,
        start_date: str = "2017-01-01",
        end_date: str = "2025-09-06",
    ) -> None:
        """Initialise the loader.

        Args:
            tickers: Yahoo Finance symbols. Defaults to BTC/ETH/SOL.
            start_date: Start of the data range.
            end_date: End of the data range.
        """
        self.tickers = tickers or list(self.TICKER_MAP.keys())
        self.start_date = start_date
        self.end_date = end_date

    def download_prices(self) -> pd.DataFrame:
        """Download daily adjusted closing prices from Yahoo Finance.

        Returns:
            DataFrame with DatetimeIndex and columns BTC, ETH, SOL.

        Raises:
            ValueError: If the downloaded DataFrame is empty.
        """
        logger.info(
            "Downloading %s  |  %s → %s",
            self.tickers,
            self.start_date,
            self.end_date,
        )
        raw = yf.download(
            self.tickers,
            start=self.start_date,
            end=self.end_date,
            auto_adjust=True,
            progress=False,
        )
        # yfinance returns MultiIndex columns when more than one ticker
        prices = (
            raw["Close"]
            if isinstance(raw.columns, pd.MultiIndex)
            else raw[["Close"]]
        )
        prices = prices.rename(columns=self.TICKER_MAP)

        if prices.empty:
            raise ValueError(
                "No price data returned. Check tickers and date range."
            )
        logger.info("Downloaded %d rows of price data.", len(prices))
        return prices

    def compute_log_returns(self, prices: pd.DataFrame) -> pd.DataFrame:
        """Compute log-returns and drop all rows containing NaN.

        Implements R: na.omit(diff(log(px))).
        Uses log-returns (not simple returns) because they are
        time-additive and better suited for VAR modelling.

        Args:
            prices: DataFrame of adjusted closing prices.

        Returns:
            DataFrame of log-returns with no missing values.
        """
        log_ret = np.log(prices).diff().dropna()
        logger.info(
            "Log-returns: %s → %s  (%d observations)",
            log_ret.index[0].date(),
            log_ret.index[-1].date(),
            len(log_ret),
        )
        return log_ret

    def load(self) -> pd.DataFrame:
        """Run full pipeline: download prices, compute log-returns.

        Returns:
            Clean log-return DataFrame ready for analysis.
        """
        prices = self.download_prices()
        return self.compute_log_returns(prices)