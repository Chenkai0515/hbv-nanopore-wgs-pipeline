# HBV Nanopore WGS Pipeline - Makefile
# Common development and deployment tasks

.PHONY: help env-dorado env-qc env-kraken2 env-daily env-medaka env-clair3 env-variants env-all
.PHONY: lint lint-check test test-fast test-all test-example clean clean-results

SHELL := /bin/bash
PYTHON := python3

# Default target
.DEFAULT_GOAL := help

help:
	@echo "═══════════════════════════════════════════════════════════════"
	@echo "              HBV Nanopore WGS Pipeline - Makefile              "
	@echo "═══════════════════════════════════════════════════════════════"
	@echo ""
	@echo "Environment Setup:"
	@echo "  env-all       Create all conda environments"
	@echo "  env-dorado    Create dorado environment (adapter trimming)"
	@echo "  env-qc        Create QC environment (FastQC, NanoPlot)"
	@echo "  env-kraken2   Create kraken2 environment (host deconv)"
	@echo "  env-daily     Create daily environment (coverage analysis)"
	@echo "  env-medaka    Create medaka environment (consensus)"
	@echo "  env-clair3    Create clair3 environment (variant calling)"
	@echo "  env-variants  Create variants environment (iVar, analysis)"
	@echo ""
	@echo "Code Quality:"
	@echo "  lint          Format Python code with black + ruff fix"
	@echo "  lint-check    Check formatting without modifying (for CI)"
	@echo ""
	@echo "Testing:"
	@echo "  test          Run fast unit tests (coordinates, utils)"
	@echo "  test-fast     Same as test (alias)"
	@echo "  test-all      Run all tests including integration"
	@echo "  test-example  Dry-run pipeline with mini dataset"
	@echo ""
	@echo "Cleanup:"
	@echo "  clean         Remove temporary files (__pycache__, .pyc, .DS_Store)"
	@echo "  clean-results Remove results directory contents"
	@echo ""
	@echo "Quick Start:"
	@echo "  make env-qc env-medaka   # Create required environments"
	@echo "  make test                # Run unit tests"
	@echo "  make test-example        # Test pipeline (dry-run)"
	@echo ""

# ═══════════════════════════════════════════════════════════════
# Environment Creation
# ═══════════════════════════════════════════════════════════════

env-dorado:
	@echo "Creating dorado environment..."
	@mamba env create -f environment/env_dorado.yml || conda env create -f environment/env_dorado.yml

env-qc:
	@echo "Creating QC environment..."
	@mamba env create -f environment/env_qc.yml || conda env create -f environment/env_qc.yml

env-kraken2:
	@echo "Creating kraken2 environment..."
	@mamba env create -f environment/env_kraken2.yml || conda env create -f environment/env_kraken2.yml

env-daily:
	@echo "Creating daily environment..."
	@mamba env create -f environment/env_daily.yml || conda env create -f environment/env_daily.yml

env-medaka:
	@echo "Creating medaka environment..."
	@mamba env create -f environment/env_medaka.yml || conda env create -f environment/env_medaka.yml

env-clair3:
	@echo "Creating clair3 environment..."
	@mamba env create -f environment/env_clair3.yml || conda env create -f environment/env_clair3.yml

env-variants:
	@echo "Creating variants environment..."
	@mamba env create -f environment/env_variants.yml || conda env create -f environment/env_variants.yml

env-all: env-dorado env-qc env-kraken2 env-daily env-medaka env-clair3 env-variants
	@echo ""
	@echo "✓ All environments created successfully!"
	@echo ""
	@echo "Activate with:"
	@echo "  conda activate env_qc"
	@echo "  conda activate env_medaka"
	@echo "  etc."

# ═══════════════════════════════════════════════════════════════
# Code Quality
# ═══════════════════════════════════════════════════════════════

lint:
	@echo "Formatting Python code..."
	@if command -v black >/dev/null 2>&1; then \
		black scripts/ tests/ --line-length 100; \
	else \
		echo "⚠ black not installed. Install with: pip install black"; \
	fi
	@if command -v ruff >/dev/null 2>&1; then \
		ruff check scripts/ tests/ --fix --ignore E501; \
	else \
		echo "⚠ ruff not installed. Install with: pip install ruff"; \
	fi
	@echo "✓ Linting complete"

lint-check:
	@echo "Checking code formatting..."
	@if command -v black >/dev/null 2>&1; then \
		black --check --diff scripts/ tests/ --line-length 100 || (echo "❌ Code needs formatting. Run 'make lint'"; exit 1); \
	else \
		echo "⚠ black not installed"; \
	fi
	@if command -v ruff >/dev/null 2>&1; then \
		ruff check scripts/ tests/ --ignore E501; \
	else \
		echo "⚠ ruff not installed"; \
	fi
	@echo "✓ Format check passed"

# ═══════════════════════════════════════════════════════════════
# Testing
# ═══════════════════════════════════════════════════════════════

# Fast unit tests (no external dependencies)
test: test-fast

test-fast:
	@echo "Running fast unit tests..."
	@if command -v pytest >/dev/null 2>&1; then \
		PYTHONPATH=scripts pytest tests/test_coordinate_transform.py tests/test_paf_utils.py tests/test_variants_utils.py -v --tb=short; \
	else \
		echo "❌ pytest not installed. Install with: pip install pytest"; \
		exit 1; \
	fi

# All tests including those that may need more setup
test-all:
	@echo "Running all tests..."
	@if command -v pytest >/dev/null 2>&1; then \
		PYTHONPATH=scripts pytest tests/ -v --tb=short; \
	else \
		echo "❌ pytest not installed. Install with: pip install pytest"; \
		exit 1; \
	fi

# Test pipeline configuration and dry-run
test-example:
	@echo "Testing pipeline with example config (dry-run)..."
	@if [ -f examples/example_config.yaml ]; then \
		echo ""; \
		echo "1. Validating config parser..."; \
		$(PYTHON) scripts/utils/config_parser.py examples/example_config.yaml --summary; \
		echo ""; \
		echo "2. Listing pipeline steps..."; \
		./run_pipeline.sh examples/example_config.yaml --list; \
		echo ""; \
		echo "3. Dry-run for sample 10090..."; \
		./run_pipeline.sh examples/example_config.yaml --sample 10090 --dry-run; \
		echo ""; \
		echo "✓ Test-example completed"; \
	else \
		echo "❌ Example config not found: examples/example_config.yaml"; \
		exit 1; \
	fi

# ═══════════════════════════════════════════════════════════════
# Cleanup
# ═══════════════════════════════════════════════════════════════

clean:
	@echo "Cleaning temporary files..."
	@find . -type f -name "*.pyc" -delete 2>/dev/null || true
	@find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	@find . -type f -name ".DS_Store" -delete 2>/dev/null || true
	@find . -type f -name "*.log" -path "*/logs/*" -delete 2>/dev/null || true
	@echo "✓ Cleaned temporary files"

clean-results:
	@echo "Cleaning results directory..."
	@if [ -d results ]; then \
		rm -rf results/*/; \
		echo "✓ Results directory cleaned"; \
	else \
		echo "No results directory found"; \
	fi
