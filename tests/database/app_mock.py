from flask import Flask


def create_app(config_name="testing"):
    """Mock app creation function for testing."""
    app = Flask(__name__)

    # Configure the app for testing
    app.config.update(
        TESTING=True,
        DUCKDB_PATH="test_db.duckdb",
        EXAMPLE_FASTA="test.faa",
        REPSEQ_FASTA="test_repseq.faa",
    )

    # Add a simple route for testing
    @app.route("/load_example", methods=["GET"])
    def load_example():
        from flask import jsonify

        try:
            with open(app.config["EXAMPLE_FASTA"], "r") as f:
                example_sequence = f.read().strip()
            return jsonify({"sequence": example_sequence})
        except Exception as e:
            app.logger.error(f"Error loading example: {str(e)}")
            return jsonify({"error": "Failed to load example sequence"}), 500

    return app
