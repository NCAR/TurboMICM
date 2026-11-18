"""Run a copied generator in isolated mode with only its own folder on sys.path.

    python3 -I run_copy.py TOOLS_DIR SCRIPT [ARGS...]
"""
import runpy, sys
from pathlib import Path

tools = Path(sys.argv[1]).resolve()
script = tools / sys.argv[2]
sys.path.insert(0, str(tools))
sys.argv = [str(script)] + sys.argv[3:]
runpy.run_path(str(script), run_name="__main__")
