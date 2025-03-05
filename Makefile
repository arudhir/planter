.PHONY: \
	help \
	up image tag push bash \
	test test-all test-verbose test-slow coverage \
	format lint check reformat \
	setup install-hooks \
	dev-setup \
	deploy \
	clean clean-build clean-pyc \
	dist release

# --- Environment---

# use Bash with brace expansion
.SHELLFLAGS = -cB
SHELL = /bin/bash

PROJECT_SLUG = planter
BUILD_IMAGE ?= $(PROJECT_SLUG)

CI_PROJECT_PATH ?= ngs-analysis/$(PROJECT_SLUG)
APP_HOME ?= /usr/src/$(PROJECT_SLUG)

help:
	@echo "Available targets:"
	@echo "  setup        - Create virtual environment and install dependencies"
	@echo "  clean        - Remove temporary files"
	@echo "  lint         - Run linters"
	@echo "  format       - Format code with Black and isort"
	@echo "  check        - Check code style without modifying files"
	@echo "  reformat     - Run format, check, and lint in sequence"
	@echo
	@echo "Testing:"
	@echo "  test         - Run fast tests"
	@echo "  test-verbose - Run tests with verbose output"
	@echo "  test-slow    - Run only slow tests"
	@echo "  test-all     - Run all tests"
	@echo "  coverage     - Run tests with coverage report"
	@echo
	@echo "Development:"
	@echo "  install      - Install editable package"
	@echo "  dev-setup    - Set up development environment with additional tools"
	@echo "  install-hooks - Install git pre-commit hooks"
	@echo
	@echo "Docker:"
	@echo "  image        - Build the image"
	@echo "  up           - Start the service"
	@echo "  tag          - Tag the image"
	@echo "  push         - Push the image to the docker registry"
	@echo "  bash         - Start a bash shell inside the docker container"
	@echo
	@echo "Packaging:"
	@echo "  dist         - Create distribution"
	@echo "  release      - Release with twine"
	@echo
	@echo "Common workflows:"
	@echo "  make check test     - Validate code style and run tests"
	@echo "  make reformat       - Format code and run linters"
	#
# --- External execution ---

MAKE_EXT = docker-compose run --rm ${PROJECT_SLUG} make -C $(APP_HOME)

# Generically execute make targets from outside the Docker container
%-ext:
	$(MAKE_EXT) $*

# --- Python ---

venv:
	source venv/bin/activate

dist: venv clean
	python -m build

install: dist
	pip install -e .

release: dist
	twine upload dist/* --verbose

# --- Cleanup ---

clean: clean-build clean-pyc

clean-build:
	rm -fr build/
	rm -fr dist/
	rm -fr *.egg-info

clean-pyc:
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +

# --- Docker image ---

image:
	docker-compose build
	#docker build --pull -t $(BUILD_IMAGE) .
	@echo built planter image as $(BUILD_IMAGE)

tag: image
	docker tag $(BUILD_IMAGE) $(TAG_IMAGE)
	@echo tagged planter image as $(TAG_IMAGE)

push: push_image = $(if $(TAG_IMAGE),$(TAG_IMAGE),$(BUILD_IMAGE))
push:
	docker push $(push_image)
	@echo pushed planter image as $(push_image)

deploy: image
	docker tag \
		$(BUILD_IMAGE) \
		docker.arudhir.com/$(CI_PROJECT_PATH)
	@echo tagged planter image
	docker push \
		docker.arudhir.com/$(CI_PROJECT_PATH)
	@echo pushed planter image
	
# --- Setup ---

setup:
	python -m pip install --upgrade pip setuptools wheel
	python -m pip install -e .
	python -m pip install pytest pytest-cov black isort flake8 mypy pre-commit

dev-setup: setup
	python -m pip install ipython jupyter notebook build twine

install-hooks:
	pre-commit install

# --- Code Quality ---

lint:
	flake8 planter
	mypy planter

format:
	black planter tests
	isort planter tests

check:
	black --check planter tests
	isort --check planter tests

reformat: format lint check

# --- Testing ---

test:
	python -m pytest tests/ -v -k "not slow"

test-verbose:
	python -m pytest tests/ -vv -s -k "not slow"

test-slow:
	python -m pytest tests/ -v -k "slow"

test-all:
	python -m pytest tests/ -v

coverage:
	python -m pytest --cov=planter --cov-report=html --cov-report=term tests/ -k "not slow"
	#
# --- Running ---

bash:
	docker-compose run --rm planter bash

server:
	FLASK_APP=app/main.py flask run --host=0.0.0.0 --port=8888
#vim: set noet
