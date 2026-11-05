#!/usr/bin/env python3
"""Build the static guide, executable-derived CLI reference, and rustdoc.

Requires Python 3, Pandoc and Cargo. Output and build caches must be outside
the source checkout. No deployment occurs. Temporary artifacts are removed
after the build; the chosen output directory must be new or empty.
"""

import argparse
import html
from html.parser import HTMLParser
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from urllib.parse import unquote, urlsplit, urlunsplit

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY = "https://github.com/PerAlexandersson/polytool/blob/main/"


def run(*args, **kwargs):
    return subprocess.run(args, check=True, text=True, **kwargs)


def captured(*args):
    return run(*args, stdout=subprocess.PIPE).stdout


def rewrite_links(node, source, destination, pages, output):
    if isinstance(node, dict):
        if node.get("t") in ("Link", "Image"):
            target = node["c"][-1][0]
            url = urlsplit(target)
            if not url.scheme and not url.netloc and url.path:
                local = (source.parent / unquote(url.path)).resolve()
                if local in pages:
                    path = os.path.relpath(output / pages[local], destination.parent)
                    node["c"][-1][0] = urlunsplit(("", "", path, url.query, url.fragment))
                elif local.is_relative_to(ROOT) and local.exists():
                    node["c"][-1][0] = REPOSITORY + local.relative_to(ROOT).as_posix()
                    if url.fragment:
                        node["c"][-1][0] += "#" + url.fragment
                else:
                    raise ValueError(f"unresolved documentation link in {source.name}: {target}")
        for value in node.values():
            rewrite_links(value, source, destination, pages, output)
    elif isinstance(node, list):
        for value in node:
            rewrite_links(value, source, destination, pages, output)


def verify_links(output):
    class Links(HTMLParser):
        def handle_starttag(self, tag, attrs):
            for key, value in attrs:
                if key in ("href", "src") and value:
                    url = urlsplit(value)
                    if not url.scheme and not url.netloc and url.path:
                        target = (self.page.parent / unquote(url.path)).resolve()
                        if not target.is_relative_to(output) or not target.exists():
                            raise ValueError(f"broken local link: {self.page}: {value}")

    parser = Links()
    # Rustdoc owns its internal links; check the guide's links into it too.
    for page in [output / "index.html", output / "README.html", *output.glob("doc/*.html")]:
        parser.page = page
        parser.feed(page.read_text())
    if any(p.is_symlink() for p in output.rglob("*")):
        raise ValueError("site artifacts must not contain symbolic links")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    workspace = Path(json.loads(captured(
        "cargo", "metadata", "--locked", "--no-deps", "--format-version", "1",
        "--manifest-path", str(ROOT / "Cargo.toml"))) ["workspace_root"])
    if output.is_relative_to(workspace) or workspace.is_relative_to(output):
        parser.error("choose an output directory outside the source workspace")
    if output.exists() and any(output.iterdir()):
        parser.error("output directory must be new or empty")
    for executable in ("cargo", "pandoc"):
        if not shutil.which(executable):
            parser.error(f"missing build tool: {executable}")
    output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="polytool-docs-", dir="/tmp") as temporary:
        temp = Path(temporary)
        env = dict(os.environ, CARGO_TARGET_DIR=str(temp / "target"), RUSTDOCFLAGS="-D warnings")
        manifest = str(ROOT / "Cargo.toml")
        run("cargo", "build", "--locked", "--manifest-path", manifest,
            "--bin", "polytool", env=env)
        # Fresh target directory prevents unrelated workspace API pages leaking
        # into the site. --no-deps retains only these three public libraries.
        run("cargo", "doc", "--locked", "--no-deps", "--lib",
            "--manifest-path", manifest, "-p", "polytool", "-p", "exact-real-roots",
            "-p", "ore-algebra", env=env)
        shutil.copytree(temp / "target/doc", output / "api")
        binary = str(temp / "target/debug/polytool")
        metadata = json.loads(captured(binary, "help", "--json"))
        reference = "# Command reference\n\nGenerated from `polytool --help` at this build.\n\n"
        reference += "```text\n" + captured(binary, "--help") + "```\n\n"
        for command in metadata["commands"]:
            name = command["name"]
            reference += f"## {name}\n\n```text\n{captured(binary, name, '--help')}```\n\n"
        for name in ("add", "check", "evidence"):
            reference += f"## oeis {name}\n\n```text\n{captured(binary, 'oeis', name, '--help')}```\n\n"
        sources = {p: p.read_text() for p in sorted((ROOT / "doc").glob("*.md"))}
        sources[ROOT / "README.md"] = (ROOT / "README.md").read_text()
        sources[ROOT / "doc/command-reference.md"] = reference
        pages = {p: p.relative_to(ROOT).with_suffix(".html") for p in sources}
        shutil.copyfile(ROOT / "scripts/docs.css", output / "style.css")
        for source, markdown in sources.items():
            destination = output / pages[source]
            destination.parent.mkdir(parents=True, exist_ok=True)
            ast = json.loads(run("pandoc", "--from=gfm", "--to=json", input=markdown,
                                 stdout=subprocess.PIPE).stdout)
            rewrite_links(ast, source, destination, pages, output)
            prefix = os.path.relpath(output, destination.parent)
            nav = temp / "nav.html"
            nav.write_text(f'<nav aria-label="Documentation"><a href="{prefix}/index.html">Polytool</a>'
                           f'<a href="{prefix}/doc/workflows.html">Workflows</a>'
                           f'<a href="{prefix}/doc/command-reference.html">Commands</a>'
                           f'<a href="{prefix}/api/polytool/index.html">Rust API</a></nav>')
            title = markdown.splitlines()[0].lstrip("# ")
            run("pandoc", "--from=json", "--to=html5", "--standalone", "--toc",
                "--metadata", f"pagetitle={title} — Polytool", "--css", f"{prefix}/style.css",
                "--include-before-body", str(nav), "--output", str(destination),
                input=json.dumps(ast))
        version = captured(binary, "--version").strip()
        (output / "index.html").write_text(
            '<!doctype html><html lang="en"><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width,initial-scale=1">'
            '<title>Polytool documentation</title><link rel="stylesheet" href="style.css">'
            '<main><h1>Polytool</h1><p>Exact polynomial tools for combinatorial research.</p>'
            '<ul><li><a href="doc/workflows.html">Start with a CLI workflow</a></li>'
            '<li><a href="doc/command-reference.html">Browse every command and option</a></li>'
            '<li><a href="api/polytool/index.html">Search the Rust API</a></li>'
            '<li><a href="doc/library.html">Use the library</a></li>'
            '<li><a href="doc/algorithms.html">Algorithms and exactness</a></li>'
            '<li><a href="doc/mcp.html">MCP tools</a></li>'
            '<li><a href="doc/web.html">Browser interface</a></li>'
            '<li><a href="README.html">Installation and project overview</a></li></ul>'
            f'<footer>Built from {html.escape(version)}. '
            '<a href="https://github.com/PerAlexandersson/polytool">Source on GitHub</a>'
            '</footer></main></html>')
        (output / ".nojekyll").write_text("")
        (output / "commands.json").write_text(json.dumps(metadata, indent=2) + "\n")
        verify_links(output)
    print(f"Documentation built and local links checked: {output / 'index.html'}")


if __name__ == "__main__":
    main()
