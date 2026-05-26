# Contributing to this Project

First off, thank you for taking the time to contribute! Contributions from the community help keep this project robust, efficient, and useful for everyone.

By participating in this project, you agree to abide by our [Code of Conduct](CODE_OF_CONDUCT.md).

## How Can I Contribute?

### Reporting Bugs
If you find a bug or unexpected behavior in the code, please search the existing Issues to see if it has already been reported. If it hasn't, open a new issue and include:
* **A clear description** of the bug and the expected behavior.
* **A minimal reproducible example** (reprex) including sample data or inputs if applicable.
* **Environment details** (e.g., R/Python version, operating system, and package versions).

### Suggesting Enhancements
We welcome ideas for new features, optimizations, or documentation improvements! Please open an issue to discuss your ideas before writing major chunks of code, so we can ensure the enhancement aligns with the project's core goals.

### Submitting Pull Requests
Ready to contribute code or documentation? Please follow this workflow:

1. **Fork the repository** and clone it locally.
2. **Create a new branch** off the main branch for your changes (`git checkout -b feature/amazing-feature` or `git checkout -b fix/bug-description`).
3. **Implement your changes** following the code style guidelines below.
4. **Test your code** thoroughly to ensure existing functionality isn't broken.
5. **Commit your changes** with clear, descriptive commit messages.
6. **Push to your fork** and open a Pull Request against our main branch.

## Code Style Guidelines

To keep the codebase uniform, readable, and easy to maintain, all submissions must strictly adhere to the following language-specific formatters before a Pull Request can be merged.

### Python Guidelines
* **Formatter:** All Python code must be formatted using **Black**. 
* Run `black your_script.py` on your files before committing. Code that triggers formatting changes during CI check pipelines will not be merged.
* **Structure:** Keep scripts modular. Break complex data processing tasks into well-documented functions rather than executing sprawling top-level code blocks.

### R Guidelines
* **Formatter:** All R code must be formatted using **Air** (the Rust-based tidyverse formatter from Posit).
* If you use **Positron**, Air is built-in; please enable **"Format on Save"** in your workspace settings. If you are using the CLI, run the following command on your files before committing:
  ```bash
  air format .
