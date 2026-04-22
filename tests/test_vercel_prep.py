from __future__ import annotations

import importlib
import os
import tempfile
import unittest
from unittest.mock import patch

from app import create_app


class VercelPrepTests(unittest.TestCase):
    def test_public_asset_route_serves_stylesheet(self) -> None:
        app = create_app()
        client = app.test_client()

        response = client.get("/assets/styles.css")

        self.assertEqual(response.status_code, 200)
        self.assertIn("text/css", response.content_type)
        self.assertIn(b".page-shell", response.data)
        response.close()

    def test_vercel_defaults_use_tmp_cache_and_disable_local_predictors(self) -> None:
        import config as config_module

        with patch.dict(os.environ, {"VERCEL": "1"}, clear=True):
            config_module = importlib.reload(config_module)
            self.assertEqual(
                config_module.Config.ANALYSIS_WORKDIR,
                os.path.join(tempfile.gettempdir(), "protein-mutation-workbench"),
            )
            self.assertEqual(config_module.Config.ENABLED_PREDICTOR_KEYS, [])

        importlib.reload(config_module)


if __name__ == "__main__":
    unittest.main()
