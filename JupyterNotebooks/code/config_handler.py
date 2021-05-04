import json
import os


class ConfigHandler:

    def __init__(self, parent=None):
        self.parent = parent

    def load(self):
        config_file_name = os.path.join(os.path.dirname(__file__), "config.json")
        with open(config_file_name) as f:
            config = json.load(f)
        self.parent.config = config
