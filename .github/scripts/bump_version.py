#!/usr/bin/env python3
"""Bump the Version field of the DESCRIPTION file from commit message prefixes.

Commit messages are classified by their prefix:

    minor:     -> patch bump (0.20.7 -> 0.20.8)
    moderate:  -> minor bump (0.20.7 -> 0.21.0)
    major:     -> major bump (0.20.7 -> 1.0.0)

The highest level found across the commits of a push is applied.  If no commit
is classified, the DESCRIPTION file is left unchanged.
"""

import os
import re
import subprocess
import sys
from pathlib import Path

PATCH, MINOR, MAJOR = 1, 2, 3

VERSION_PATTERN = re.compile(r"(?m)^Version: (\d+)\.(\d+)\.(\d+)[ \t]*$")


def classify(message):
    """Return the bump level requested by a single commit message."""
    if re.search(r"(?im)^[ \t]*major:", message):
        return MAJOR
    if re.search(r"(?im)^[ \t]*moderate:", message):
        return MINOR
    if re.search(r"(?im)^[ \t]*minor:", message):
        return PATCH
    return 0


def highest_level(messages):
    """Return the highest bump level requested by any of the messages."""
    return max((classify(message) for message in messages), default=0)


def bump(version, level):
    """Return `version` bumped by `level`."""
    major, minor, patch = (int(part) for part in version.split("."))
    if level == MAJOR:
        return f"{major + 1}.0.0"
    if level == MINOR:
        return f"{major}.{minor + 1}.0"
    if level == PATCH:
        return f"{major}.{minor}.{patch + 1}"
    return version


def commit_messages(before, after, cwd=None):
    """Return the messages of the commits introduced by a push."""
    revision = after if not before or set(before) == {"0"} else f"{before}..{after}"
    output = subprocess.check_output(
        ["git", "log", revision, "--format=%B%x00"], text=True, cwd=cwd
    )
    return [message for message in output.split("\0") if message.strip()]


def update_description(path, level):
    """Apply `level` to the Version field of the DESCRIPTION file at `path`.

    Returns the new version, or None when nothing was changed.
    """
    if not level:
        return None
    path = Path(path)
    text = path.read_text()
    match = VERSION_PATTERN.search(text)
    if not match:
        raise SystemExit("Expected an X.Y.Z Version field in DESCRIPTION")
    version = bump(".".join(match.groups()), level)
    path.write_text(f"{text[:match.start()]}Version: {version}{text[match.end():]}")
    return version


def main():
    before = os.environ.get("BEFORE", "")
    after = os.environ.get("AFTER", "HEAD")
    description = os.environ.get("DESCRIPTION_PATH", "DESCRIPTION")

    level = highest_level(commit_messages(before, after))
    version = update_description(description, level)
    if version is None:
        print("No classified commits found; leaving the version unchanged.")
    else:
        print(f"Bumped version to {version}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
