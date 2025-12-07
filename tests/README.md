# Tests

This directory contains unit tests and integration tests for the pipeline.

## Test Structure

```
tests/
├── README.md                      # This file
├── conftest.py                    # Pytest fixtures
├── test_coordinate_transform.py   # Coordinate transformation tests
├── test_mask_hbv_regions.py       # HBV masking tests
└── ...                            # Additional tests
```

## Running Tests

### Prerequisites

```bash
# Install test dependencies
conda activate hbv_base
pip install pytest pytest-cov
```

### Run All Tests

```bash
# From repository root
pytest tests/ -v
```

### Run Specific Tests

```bash
# Single test file
pytest tests/test_coordinate_transform.py -v

# Specific test function
pytest tests/test_coordinate_transform.py::test_rotated_to_original -v
```

### With Coverage

```bash
pytest tests/ --cov=scripts --cov-report=html
# View report in htmlcov/index.html
```

## Test Categories

### Unit Tests

Test individual functions in isolation:
- Coordinate transformation formulas
- Quality score calculations
- Filter logic

### Integration Tests

Test module interactions:
- QC pipeline end-to-end
- Variant calling workflow
- Coordinate transformation chain

## Writing Tests

### Template

```python
# tests/test_example.py
import pytest

def test_example_function():
    """Test description."""
    # Arrange
    input_data = ...
    expected = ...
    
    # Act
    result = function_under_test(input_data)
    
    # Assert
    assert result == expected
```

### Fixtures

Common fixtures are in `conftest.py`:

```python
@pytest.fixture
def sample_fastq():
    """Provide sample FASTQ content for testing."""
    return "@read1\nATCGATCG\n+\nIIIIIIII\n"

@pytest.fixture
def reference_sequence():
    """Provide reference sequence."""
    return "ATCGATCGATCG..."
```

## Continuous Integration

Tests are run automatically via GitHub Actions on:
- Push to main or develop branch
- Pull requests to main

See `.github/workflows/ci.yml` for CI configuration.

Note: CI only runs pure Python tests (coordinate transforms, config parsing).
Heavy integration tests requiring bioinformatics tools should be run locally.

## Test Data

Test data is stored in:
- `tests/data/`: Small test files
- `examples/mini_dataset/`: Larger test datasets

## Adding New Tests

1. Create test file: `tests/test_module_name.py`
2. Add fixtures to `conftest.py` if needed
3. Run tests locally: `pytest tests/test_module_name.py -v`
4. Ensure all tests pass before committing

