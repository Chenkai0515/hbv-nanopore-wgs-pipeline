# HBV Nanopore WGS Pipeline - Makefile
# Common development and deployment tasks

.PHONY: help env env-all lint test test-example clean

SHELL := /bin/bash

help:
	@echo "HBV Nanopore Pipeline - Available targets:"
	@echo ""
	@echo "  env-dorado    Create dorado environment (adapter trimming)"
	@echo "  env-qc        Create QC environment (FastQC, NanoPlot)"
	@echo "  env-kraken2   Create kraken2 environment (host deconv)"
	@echo "  env-daily     Create daily environment (coverage analysis)"
	@echo "  env-medaka    Create medaka environment (consensus)"
	@echo "  env-clair3    Create clair3 environment (variant calling)"
	@echo "  env-variants  Create variants environment (iVar, analysis)"
	@echo "  env-all       Create all environments"
	@echo ""
	@echo "  lint          Format Python code with black"
	@echo "  test          Run unit tests"
	@echo "  test-example  Run pipeline on example data (dry-run)"
	@echo "  clean         Remove temporary files"
	@echo ""

# Individual environment targets
env-dorado:
	mamba env create -f environment/env_dorado.yml || conda env create -f environment/env_dorado.yml

env-qc:
	mamba env create -f environment/env_qc.yml || conda env create -f environment/env_qc.yml

env-kraken2:
	mamba env create -f environment/env_kraken2.yml || conda env create -f environment/env_kraken2.yml

env-daily:
	mamba env create -f environment/env_daily.yml || conda env create -f environment/env_daily.yml

env-medaka:
	mamba env create -f environment/env_medaka.yml || conda env create -f environment/env_medaka.yml

env-clair3:
	mamba env create -f environment/env_clair3.yml || conda env create -f environment/env_clair3.yml

env-variants:
	mamba env create -f environment/env_variants.yml || conda env create -f environment/env_variants.yml

# Create all environments
env-all: env-dorado env-qc env-kraken2 env-daily env-medaka env-clair3 env-variants
	@echo "All environments created successfully"

# Code formatting
lint:
	@if command -v black >/dev/null 2>&1; then \
		black scripts/**/*.py; \
	else \
		echo "black not installed. Install with: pip install black"; \
	fi

# Run tests
test:
	@if command -v pytest >/dev/null 2>&1; then \
		pytest tests/ -v; \
	else \
		echo "pytest not installed. Install with: pip install pytest"; \
	fi

# Test pipeline with example config (dry-run)
test-example:
	@if [ -f examples/example_config.yaml ]; then \
		./run_pipeline.sh examples/example_config.yaml --dry-run --list; \
	else \
		echo "Example config not found"; \
	fi

# Clean temporary files
clean:
	find . -type f -name "*.pyc" -delete
	find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name ".DS_Store" -delete
	@echo "Cleaned temporary files"

