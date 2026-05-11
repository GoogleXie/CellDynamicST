# Contributing to CellDynamicST

Thank you for your interest in contributing to CellDynamicST. This document provides guidelines and information for contributors.

## How to Contribute

### Reporting Bugs

If you find a bug, please open an issue on GitHub using the **Bug Report** template. Include a minimal reproducible example, the full error message, and your `sessionInfo()` output.

### Suggesting Features

Feature requests are welcome. Please open an issue using the **Feature Request** template and describe the use case, expected behavior, and any alternatives you have considered.

### Pull Requests

1. Fork the repository and create a new branch from `main`.
2. Make your changes, following the code style guidelines below.
3. Add or update tests in `tests/testthat/` for any new functionality.
4. Ensure `R CMD check` passes with no errors or warnings.
5. Update documentation (roxygen2 comments) for any changed functions.
6. Submit a pull request with a clear description of the changes.

## Code Style

CellDynamicST follows the tidyverse style guide with these conventions:

- Function names use `cdst_verb_noun()` pattern for exported functions.
- Internal functions use `snake_case` without the `cdst_` prefix.
- All exported functions must have complete roxygen2 documentation.
- Use `checkmate` for input validation in all public functions.
- Use `cli` for user-facing messages (not `message()` or `cat()`).

## Development Setup

```r
# Install development dependencies
install.packages(c("devtools", "testthat", "roxygen2", "lintr"))

# Clone and install in development mode
devtools::load_all(".")

# Run tests
devtools::test()

# Check package
devtools::check()
```

## Code of Conduct

This project follows the [Contributor Covenant Code of Conduct](CODE_OF_CONDUCT.md). By participating, you are expected to uphold this code.
