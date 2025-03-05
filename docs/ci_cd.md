# Continuous Integration and Deployment

This document describes the CI/CD pipelines for the Planter project.

## GitHub Actions Workflows

Planter uses GitHub Actions for Continuous Integration and Deployment. There are two main workflows:

### CI Workflow

The CI workflow (`ci.yml`) runs on all pushes to the main branch and on all pull requests. It consists of the following jobs:

1. **Test** - Runs all fast tests (excluding slow tests) on Python 3.12:
   - Runs the tests with pytest
   - Generates and uploads code coverage reports

2. **Lint** - Runs code quality checks:
   - Black for code formatting
   - Flake8 for lint checks
   - isort for import ordering
   - mypy for type checking

3. **Slow Tests** - Runs tests marked as slow, but only on the main branch
   - Includes end-to-end workflow tests
   - Uses Python 3.12

4. **Build Image** - Builds the Docker image for verification
   - Runs only on the main branch
   - Depends on the test and lint jobs

### Release Workflow

The release workflow (`release.yml`) runs when a new release is created on GitHub. It performs the following tasks:

1. **Build and Publish Python Package**
   - Builds the Python package (sdist and wheel)
   - Uploads to PyPI

2. **Build and Push Docker Image**
   - Builds the Docker image
   - Tags with release version and "latest"
   - Pushes to Docker registry

## Working with the CI System

### Running Tests Locally

To run the same tests that CI will run, you can use the following commands:

```bash
# Run fast tests (what CI runs on PRs)
make test

# Run all tests (including slow tests)
make test-all

# Run only slow tests
python -m pytest tests/ -v -k "slow"
```

### Adding New Slow Tests

When adding new tests that might take a long time to run (more than a few seconds), mark them with `@pytest.mark.slow`:

```python
import pytest

@pytest.mark.slow
def test_my_slow_test():
    # Test code here
```

This ensures the test is skipped in regular CI runs and only executed in the dedicated slow test job.

### Code Coverage

The CI system automatically generates code coverage reports. These are uploaded to Codecov, which tracks code coverage over time and displays coverage reports for pull requests.

## Deploying New Releases

To deploy a new release:

1. Update version in `pyproject.toml`
2. Create a new release in GitHub with a tag that matches the version
3. The release workflow will automatically build and publish the package to PyPI and push the Docker image.