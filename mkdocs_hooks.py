"""MkDocs hooks for values sourced from the package metadata."""

import ast
from pathlib import Path


def _read_package_version():
    init_path = Path(__file__).resolve().parent / "crystmatch" / "__init__.py"
    module = ast.parse(init_path.read_text(encoding="utf-8"), filename=str(init_path))

    for node in module.body:
        if not isinstance(node, ast.Assign):
            continue
        if any(isinstance(target, ast.Name) and target.id == "__version__" for target in node.targets):
            return ast.literal_eval(node.value)

    raise RuntimeError(f"Could not find __version__ in {init_path}")


PACKAGE_VERSION = _read_package_version()


def on_page_markdown(markdown, **kwargs):
    """Replace the package-version placeholder before Markdown rendering."""
    return markdown.replace("{{ crystmatch_version }}", PACKAGE_VERSION)
