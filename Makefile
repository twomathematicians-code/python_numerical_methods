.PHONY: all run run-all list test clean help

METHOD ?= 01_bisection

# Default: show help
help:
	@echo "Python Numerical Methods Toolkit"
	@echo "==============================="
	@echo ""
	@echo "  make list             List all available methods"
	@echo "  make run METHOD=XX    Run a specific method (e.g. make run METHOD=01_bisection)"
	@echo "  make run-all          Run all 20 methods with demo inputs"
	@echo "  make test             Run the automated test suite"
	@echo "  make clean            Remove Python cache files"
	@echo ""

list:
	@echo "Available numerical methods:"
	@ls src/*.py | sed 's/src\///' | sed 's/.py//'

run:
	@python src/$(METHOD).py

run-all:
	@for f in src/*.py; do echo ""; echo ">>> Running $$(basename $$f)"; python "$$f" < /dev/null 2>&1 || true; done

test:
	python tests/test_all_methods.py

clean:
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -name "*.pyc" -delete 2>/dev/null || true
