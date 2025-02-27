.PHONY: \
	help \
	up image tag push bash \
	test test-all \
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
	@echo
	@echo "Usage: make [target]"
	@echo "n.b., add -ext to run target in container"
	@echo
	@echo "Cleanup:"
	@echo "    clean              clean-build, clean-docs, clean-node, clean-coverage, and clean-pyc"
	@echo "    clean-build        remove build artifacts"
	@echo "    clean-pyc          remove Python file artifacts"
	@echo
	@echo "Docker image:"
	@echo "    up                 start the service"
	@echo "    image              build the image"
	@echo "    tag                tag the image"
	@echo "    push               push the image to the docker registry"
	@echo "    bash               start a bash shell inside the docker container"
	@echo
	@echo "Testing:"
	@echo "    test               run fast tests"
	@echo "    test-all           run all tests"
	@echo
	@echo "Packaging:"
	@echo "    venv               create virtualenv"
	@echo "    dist               create distribution"
	@echo "    install            install editable package"
	@echo "    release            release with twine"
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
	
# --- Testing ---

test:
	python -m pytest tests/ -v

test-all:
	pytest --runslow --html=pytest-full.html --self-contained-html
	#
# --- Running ---

bash:
	docker-compose run --rm planter bash

server:
	FLASK_APP=app/main.py flask run --host=0.0.0.0 --port=8888
#vim: set noet
