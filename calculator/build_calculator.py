"""Build the single-file web calculator: insert doac_core.js into template.html -> index.html.
Usage: python calculator/build_calculator.py
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
core = open(os.path.join(HERE, "doac_core.js"), encoding="utf-8").read()
tpl = open(os.path.join(HERE, "template.html"), encoding="utf-8").read()
assert "/*CORE*/" in tpl, "placeholder /*CORE*/ not found in template.html"
open(os.path.join(HERE, "index.html"), "w", encoding="utf-8").write(tpl.replace("/*CORE*/", core))
print("calculator/index.html written")
