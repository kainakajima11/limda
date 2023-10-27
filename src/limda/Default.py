import pathlib 
import yaml
from typing import List, Dict, Final

DOT_PATH : Final[pathlib.Path] = pathlib.Path.home() / ".limda.yaml"
DEFAULT : Final[Dict] = {}

PARA : Final[List[str]] = None

def load_yaml(path : pathlib.Path = None)->Dict:
    """
    yamlfileを読み込む
    """
    with open(path, "r") as y:
        default = yaml.safe_load(y)
    return default

def set_default():
    """
    ~/.limda.yaml から default値を設定する.
    """
    if not pathlib.Path.exists(DOT_PATH):
        print("NO DEFAULT")
        return

    global DEFAULT
    DEFAULT = load_yaml(DOT_PATH)

    if "PARA" in DEFAULT:
        global PARA
        PARA = DEFAULT["PARA"] 