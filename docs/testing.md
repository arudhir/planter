# Planter Testing Guide

This guide explains how to run and write tests for the Planter project.

## Running Tests

### Basic Commands

```bash
# Run all tests
make test

# Or using pytest directly
python -m pytest tests/ -v

# Run a specific test file
python -m pytest tests/database/test_schema_version.py -v

# Run a specific test class
python -m pytest tests/database/test_schema_version.py::TestSchemaVersion -v

# Run a specific test method
python -m pytest tests/database/test_schema_version.py::TestSchemaVersion::test_get_db_schema_version -v

# Show output during tests (useful for debugging)
python -m pytest tests/ -v -s
```

### Running With Coverage

```bash
# Generate coverage report
python -m pytest --cov=planter --cov-report=term tests/

# Generate HTML coverage report
python -m pytest --cov=planter --cov-report=html tests/
```

## Test Structure

The Planter codebase follows standard Python unittest practices:

1. Test files are named with the prefix `test_`
2. Test classes inherit from `unittest.TestCase`
3. Test methods have names starting with `test_`

### Example Test Class

```python
class TestSequenceDBBuilder(unittest.TestCase):
    """Test cases for the SequenceDBBuilder class."""

    def setUp(self):
        """Set up test environment before each test."""
        self.temp_dir = tempfile.mkdtemp()
        self.db_path = os.path.join(self.temp_dir, "test.duckdb")
        # More setup...

    def tearDown(self):
        """Clean up after each test."""
        shutil.rmtree(self.temp_dir)

    def test_database_initialization(self):
        """Test that the database is properly initialized."""
        # Test code...
        self.assertEqual(len(sequences), 2)
        self.assertIn('sequences', table_names)
```

## Test Fixtures

Test fixture data is stored in `tests/fixtures/` and includes:

- Sample database files
- Test sequence data
- Example configuration files

When using fixtures, always use relative paths and the `Path` object from the standard library:

```python
from pathlib import Path

# Get path to fixtures directory
fixtures_dir = Path(__file__).parent.parent / "fixtures"
test_file = fixtures_dir / "SRR12068547" / "SRR12068547.duckdb"
```

## Mocking

For tests that involve external services or expensive operations, we use `unittest.mock`:

```python
from unittest.mock import patch, MagicMock

@patch('planter.database.builder.get_sra_info')
def test_database_initialization(self, mock_get_sra_info):
    # Set up the mock
    mock_get_sra_info.return_value = {
        'organism': 'Test Organism',
        'study_title': 'Test Study',
        # Other mock data...
    }
    
    # Test code...
```

## Writing New Tests

When writing new tests:

1. **Identify the component** to test (function, class, module)
2. **Create necessary fixtures** or use existing ones
3. **Write test methods** covering various scenarios:
   - Normal operation
   - Edge cases
   - Error conditions
4. **Use appropriate assertions** (assertEqual, assertTrue, assertRaises, etc.)
5. **Keep tests independent** - each test should be runnable on its own

## Test Categories

The Planter codebase has several test categories:

### Database Tests

Tests for database operations, schema management, and queries:
- `tests/database/test_database_builder.py`
- `tests/database/test_schema_version.py`
- `tests/database/test_expression_queries.py`

### Pipeline Tests

Tests for the Snakemake workflow pipeline:
- `tests/pipeline/test_duckdb_creation.py`

### Utility Tests

Tests for utility functions:
- `tests/database/utils/test_duckdb_utils.py`
- `tests/database/utils/test_s3_utils.py`

## Continuous Integration

Tests are automatically run in CI on pull requests and merges to main. Make sure all tests pass before submitting PRs.