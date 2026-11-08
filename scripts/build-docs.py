#!/usr/bin/env python3
"""Build the static guide, executable-derived CLI reference, and rustdoc.

Requires Python 3, Pandoc and Cargo. Output and build caches must be outside
the source checkout. No deployment occurs. Temporary artifacts are removed
after the build; the chosen output directory must be new or empty.

Pages: every doc/*.md guide (new guides such as doc/python.md are picked up
automatically), README.md, the command reference, and rustdoc under api/.
The command reference has two sources. If `polytool help --json` describes
commands richly (entries carry usage, description, options or examples, and
optionally a top-level "topics" list), one page per command and topic is
generated under doc/commands/ and doc/topics/ (relative links in their text
resolve from those directories). Otherwise, or with
`--command-source scrape`, the reference is scraped from each command's
`--help` text.
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


# Keys that mark a `help --json` command entry as rich enough to drive pages.
RICH_KEYS = ("usage", "description", "options", "examples", "exit_codes")


def slug(name):
    return "-".join("".join(c if c.isalnum() else " " for c in name.lower()).split())


def load_help_metadata(binary):
    """Return the parsed `help --json` document, or None if unsupported."""
    result = subprocess.run([binary, "help", "--json"], text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    if result.returncode != 0:
        return None
    try:
        metadata = json.loads(result.stdout)
    except json.JSONDecodeError:
        return None
    if not isinstance(metadata, dict) or not isinstance(metadata.get("commands"), list):
        return None
    return metadata


def is_rich(metadata):
    commands = metadata["commands"]
    return bool(commands) and all(
        isinstance(c, dict) and any(k in c for k in RICH_KEYS) for c in commands)


def flatten_commands(commands, prefix=""):
    """Yield (full name, entry) including nested "subcommands"."""
    for command in commands:
        name = (prefix + " " + command["name"]).strip()
        yield name, command
        yield from flatten_commands(command.get("subcommands", []), name)


def render_option(option):
    if isinstance(option, str):
        return f"- `{option}`\n"
    flags = option.get("flags") or option.get("name") or ""
    if isinstance(flags, list):
        flags = ", ".join(flags)
    value = option.get("value")
    label = f"{flags} {value}" if value else flags
    text = f"- `{label}`"
    if option.get("description"):
        text += " \u2014 " + str(option["description"]).strip()
    if option.get("default") is not None:
        text += f" (default: `{option['default']}`)"
    return text + "\n"


def command_page(name, command):
    page = f"# {name}\n\n"
    if command.get("summary"):
        page += f"{command['summary']}\n\n"
    if command.get("aliases"):
        page += "Aliases: " + ", ".join(f"`{a}`" for a in command["aliases"]) + "\n\n"
    if command.get("usage"):
        usage = command["usage"]
        usage = "\n".join(usage) if isinstance(usage, list) else usage
        page += f"## Usage\n\n```text\n{usage}\n```\n\n"
    if command.get("description"):
        page += f"{str(command['description']).strip()}\n\n"
    if command.get("options"):
        page += "## Options\n\n" + "".join(render_option(o) for o in command["options"]) + "\n"
    if command.get("examples"):
        page += "## Examples\n\n"
        for example in command["examples"]:
            text = example if isinstance(example, str) else example.get("command", "")
            note = "" if isinstance(example, str) else example.get("description", "")
            page += (note + "\n\n" if note else "") + f"```sh\n{text}\n```\n\n"
    if command.get("exit_codes"):
        page += "## Exit codes\n\n"
        codes = command["exit_codes"]
        if isinstance(codes, dict):
            page += "".join(f"- `{k}` \u2014 {v}\n" for k, v in codes.items())
        else:
            page += "".join(f"- {c}\n" for c in codes)
        page += "\n"
    return page


def rich_sources(binary, metadata):
    """Markdown sources for the command index and one page per command/topic."""
    sources = {}
    index = "# Command reference\n\nGenerated from `polytool help --json` at this build.\n\n"
    index += "```text\n" + captured(binary, "--help") + "```\n\n## Commands\n\n"
    group = None
    for name, command in flatten_commands(metadata["commands"]):
        if command.get("group") != group and command.get("group"):
            group = command["group"]
            index += f"\n### {group}\n\n"
        path = ROOT / "doc/commands" / (slug(name) + ".md")
        sources[path] = command_page(name, command)
        summary = f" \u2014 {command['summary']}" if command.get("summary") else ""
        index += f"- [`{name}`](commands/{slug(name)}.md){summary}\n"
    topics = [t for t in metadata.get("topics", []) if isinstance(t, dict) and t.get("name")]
    if topics:
        index += "\n## Topics\n\n"
        for topic in topics:
            path = ROOT / "doc/topics" / (slug(topic["name"]) + ".md")
            title = topic.get("title") or topic["name"]
            body = topic.get("body") or topic.get("description") or ""
            sources[path] = f"# {title}\n\n{str(body).strip()}\n"
            index += f"- [{title}](topics/{slug(topic['name'])}.md)\n"
    sources[ROOT / "doc/command-reference.md"] = index
    return sources


def scraped_sources(binary, metadata):
    """Command reference scraped from each command's --help text."""
    reference = "# Command reference\n\nGenerated from `polytool --help` at this build.\n\n"
    reference += "```text\n" + captured(binary, "--help") + "```\n\n"
    for command in (metadata or {}).get("commands", []):
        name = command["name"]
        reference += f"## {name}\n\n```text\n{captured(binary, name, '--help')}```\n\n"
    for name in ("add", "check", "evidence"):
        reference += f"## oeis {name}\n\n```text\n{captured(binary, 'oeis', name, '--help')}```\n\n"
    return {ROOT / "doc/command-reference.md": reference}


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
    for page in [output / "index.html", output / "README.html", *output.glob("doc/**/*.html")]:
        parser.page = page
        parser.feed(page.read_text())
    if any(p.is_symlink() for p in output.rglob("*")):
        raise ValueError("site artifacts must not contain symbolic links")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--command-source", choices=("auto", "scrape"), default="auto",
                        help="auto: use `help --json` pages when it is rich, else scrape "
                             "--help; scrape: always scrape")
    parser.add_argument("--target-dir", type=Path,
                        help="reusable Cargo target directory outside the workspace "
                             "(default: a fresh temporary one); its doc/ is rebuilt each run")
    parser.add_argument("--binary", type=Path,
                        help="use this polytool binary for help output instead of building one "
                             "(used to test the command-reference paths)")
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
        target = args.target_dir.resolve() if args.target_dir else temp / "target"
        if target.is_relative_to(workspace):
            parser.error("choose a --target-dir outside the source workspace")
        env = dict(os.environ, CARGO_TARGET_DIR=str(target), RUSTDOCFLAGS="-D warnings")
        manifest = str(ROOT / "Cargo.toml")
        if not args.binary:
            run("cargo", "build", "--locked", "--manifest-path", manifest,
                "--bin", "polytool", env=env)
        shutil.rmtree(target / "doc", ignore_errors=True)
        # Fresh target directory prevents unrelated workspace API pages leaking
        # into the site. --no-deps retains only these three public libraries.
        run("cargo", "doc", "--locked", "--no-deps", "--lib",
            "--manifest-path", manifest, "-p", "polytool", "-p", "exact-real-roots",
            "-p", "ore-algebra", env=env)
        shutil.copytree(target / "doc", output / "api")
        binary = str(args.binary.resolve() if args.binary else target / "debug/polytool")
        metadata = load_help_metadata(binary)
        rich = args.command_source == "auto" and metadata is not None and is_rich(metadata)
        print("Command reference source:",
              "help --json (per-command pages)" if rich else "scraped --help")
        sources = {p: p.read_text() for p in sorted((ROOT / "doc").glob("*.md"))}
        sources[ROOT / "README.md"] = (ROOT / "README.md").read_text()
        sources.update(rich_sources(binary, metadata) if rich else scraped_sources(binary, metadata))
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
            python = (f'<a href="{prefix}/doc/python.html">Python</a>'
                      if ROOT / "doc/python.md" in sources else "")
            nav.write_text(f'<nav aria-label="Documentation"><a href="{prefix}/index.html">Polytool</a>'
                           f'<a href="{prefix}/doc/workflows.html">Workflows</a>'
                           f'<a href="{prefix}/doc/command-reference.html">Commands</a>'
                           f'{python}'
                           f'<a href="{prefix}/api/polytool/index.html">Rust API</a></nav>')
            title = markdown.splitlines()[0].lstrip("# ")  # first line is the heading
            run("pandoc", "--from=json", "--to=html5", "--standalone", "--toc",
                "--metadata", f"pagetitle={title} — Polytool", "--css", f"{prefix}/style.css",
                "--include-before-body", str(nav), "--output", str(destination),
                input=json.dumps(ast))
        version = captured(binary, "--version").strip()
        featured = {"workflows", "command-reference", "library", "algorithms", "mcp", "web",
                    "python", "development"}
        more = "".join(
            f'<li><a href="{pages[p].as_posix()}">{html.escape(sources[p].splitlines()[0].lstrip("# "))}</a></li>'
            for p in sorted(sources, key=lambda q: pages[q].as_posix())
            if p.parent == ROOT / "doc" and p.stem not in featured)
        python_item = ('<li><a href="doc/python.html">Python bindings</a></li>'
                       if ROOT / "doc/python.md" in sources else "")
        development_item = ('<li><a href="doc/development.html">Development: tests, CI and docs</a></li>'
                            if ROOT / "doc/development.md" in sources else "")
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
            f'{python_item}{development_item}'
            '<li><a href="README.html">Installation and project overview</a></li></ul>'
            f'<h2>More guides</h2><ul>{more}</ul>'
            f'<footer>Built from {html.escape(version)}. '
            '<a href="https://github.com/PerAlexandersson/polytool">Source on GitHub</a>'
            '</footer></main></html>')
        (output / ".nojekyll").write_text("")
        if metadata is not None:
            (output / "commands.json").write_text(json.dumps(metadata, indent=2) + "\n")
        verify_links(output)
    print(f"Documentation built and local links checked: {output / 'index.html'}")


if __name__ == "__main__":
    main()
