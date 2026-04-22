from __future__ import annotations

import os
from pathlib import Path

from config import Config


def create_app():
    from flask import Flask

    _load_predictor_env_file()

    app = Flask(__name__)
    app.config.from_object(Config)

    from app.routes.main import bp as main_bp

    app.register_blueprint(main_bp)
    return app


def _load_predictor_env_file() -> None:
    env_path = Path(__file__).resolve().parents[1] / ".env.predictors"
    if not env_path.exists():
        return

    for raw_line in env_path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue

        key, value = line.split("=", 1)
        key = key.strip()
        if key and key not in os.environ:
            os.environ[key] = value.strip()
