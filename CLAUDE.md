# Planter Build and Test Commands

## Build Commands
- Build Docker image: `make image`
- Install package locally: `pip install -e .`

## Test Commands
- Run all tests: `make test` or `python -m pytest tests/ -v`
- Run a single test file: `python -m pytest tests/path/to/test_file.py -v`
- Run a specific test: `python -m pytest tests/path/to/test_file.py::TestClass::test_method -v`
- Show test output: `python -m pytest tests/ -v -s`

## Code Style Guidelines
- Use type hints for function parameters and return values
- Follow PEP 8 naming conventions: snake_case for variables/functions, CamelCase for classes
- Use dataclasses for data containers
- Include docstrings for all modules, classes, and methods
- Use context managers for resource management (with statement)
- Use pathlib.Path for file path handling instead of os.path
- Proper error handling with specific exceptions and informative messages
- Organize imports: standard library, then third-party, then local modules