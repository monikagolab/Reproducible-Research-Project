# Reproducible-Research-Project
This repository contains materials prepared as part of the Reproducible Research project for course conducted during the summer semester of the 2025/2026 academic year.

## Project Structure

```
Reproducible-Research-Project/
├── crypto_var/               # Main Python package
│   ├── __init__.py
│   ├── data.py               # Data download & log-returns (CryptoDataLoader)
│   ├── model.py              # VAR estimation & diagnostics
│   ├── forecast.py           # Out-of-sample validation & metrics
│   └── plots.py              # All visualizations
├── notebooks/
│   └── exploration.ipynb     # Development/scratch notebook
├── report/
│   └── report.qmd            # Quarto report (final output)
├── docs/                     # Sphinx HTML documentation (auto-generated)
├── tests/                    # Unit tests
├── Dockerfile                # Docker image definition
├── docker-compose.yml        # Docker Compose configuration
├── Makefile                  # Automation of key project steps
├── pyproject.toml            # Linting configuration (ruff)
├── .pre-commit-config.yaml   # Pre-commit hooks
├── requirements.txt          # Python dependencies
└── README.md
```
