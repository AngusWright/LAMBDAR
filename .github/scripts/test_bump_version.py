#!/usr/bin/env python3
"""Tests for the version bumping workflow and its helper script."""

import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import yaml

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent
SCRIPT = SCRIPT_DIR / "bump_version.py"
WORKFLOW = REPO_ROOT / ".github" / "workflows" / "version.yml"

sys.path.insert(0, str(SCRIPT_DIR))

import bump_version  # noqa: E402

DESCRIPTION_TEMPLATE = """Package: LAMBDAR
Type: Package
Version: {version}
Date: 2020-11-12
License: GPL-2
"""


class WorkflowStructureTest(unittest.TestCase):
    """Initial setup test: the workflow exists and is wired up correctly."""

    def setUp(self):
        self.workflow = yaml.safe_load(WORKFLOW.read_text())

    def test_runs_on_pushes_to_master(self):
        # `on` is parsed as the boolean True by the YAML 1.1 parser.
        triggers = self.workflow[True]
        self.assertEqual(triggers["push"]["branches"], ["master"])

    def test_requests_write_permission(self):
        self.assertEqual(self.workflow["permissions"]["contents"], "write")

    def test_runs_the_bump_script_and_pushes(self):
        steps = self.workflow["jobs"]["bump"]["steps"]
        run_steps = " ".join(step.get("run", "") for step in steps)
        self.assertIn("python3 .github/scripts/bump_version.py", run_steps)
        self.assertIn("git push", run_steps)

    def test_checks_out_the_full_history(self):
        checkout = self.workflow["jobs"]["bump"]["steps"][0]
        self.assertEqual(checkout["with"]["fetch-depth"], 0)


class BumpTest(unittest.TestCase):
    """Unit tests for the version arithmetic and commit classification."""

    def test_patch_bump(self):
        self.assertEqual(bump_version.bump("0.20.7", bump_version.PATCH), "0.20.8")

    def test_minor_bump(self):
        self.assertEqual(bump_version.bump("0.20.7", bump_version.MINOR), "0.21.0")

    def test_major_bump(self):
        self.assertEqual(bump_version.bump("0.20.7", bump_version.MAJOR), "1.0.0")

    def test_highest_level_wins(self):
        messages = ["minor: a", "major: b", "moderate: c"]
        self.assertEqual(bump_version.highest_level(messages), bump_version.MAJOR)

    def test_unclassified_messages(self):
        self.assertEqual(bump_version.highest_level(["fix a typo"]), 0)
        self.assertEqual(bump_version.highest_level([]), 0)

    def test_prefix_must_start_a_line(self):
        self.assertEqual(bump_version.highest_level(["see major: notes"]), 0)
        self.assertEqual(
            bump_version.highest_level(["summary\n\nmajor: body"]), bump_version.MAJOR
        )


class EndToEndTest(unittest.TestCase):
    """Run the script against a throwaway git repository, as the workflow does."""

    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.directory = Path(directory.name)
        self.description = self.directory / "DESCRIPTION"
        self.description.write_text(DESCRIPTION_TEMPLATE.format(version="0.20.7"))
        self.git("init", "-q", "-b", "master")
        self.git("config", "user.name", "Test")
        self.git("config", "user.email", "test@example.com")
        self.commit("initial commit")
        self.before = self.head()

    def git(self, *arguments):
        return subprocess.check_output(
            ["git", *arguments], cwd=self.directory, text=True
        )

    def commit(self, message):
        self.git("add", "-A")
        self.git("commit", "-q", "--allow-empty", "-m", message)

    def head(self):
        return self.git("rev-parse", "HEAD").strip()

    def run_script(self):
        environment = dict(
            os.environ,
            BEFORE=self.before,
            AFTER=self.head(),
            DESCRIPTION_PATH=str(self.description),
        )
        subprocess.check_call(
            [sys.executable, str(SCRIPT)], cwd=self.directory, env=environment
        )
        return self.version()

    def version(self):
        match = bump_version.VERSION_PATTERN.search(self.description.read_text())
        return match.group(0).split(": ")[1]

    def test_patch_bump(self):
        self.commit("minor: tidy up the aperture code")
        self.assertEqual(self.run_script(), "0.20.8")

    def test_minor_bump(self):
        self.commit("moderate: add a new deprojected CoG option")
        self.assertEqual(self.run_script(), "0.21.0")

    def test_major_bump(self):
        self.commit("major: change the catalogue interface")
        self.assertEqual(self.run_script(), "1.0.0")

    def test_highest_level_of_a_push_is_applied(self):
        self.commit("minor: tidy up")
        self.commit("major: breaking change")
        self.commit("moderate: new feature")
        self.assertEqual(self.run_script(), "1.0.0")

    def test_no_classified_commits_leaves_the_version_unchanged(self):
        self.commit("update the documentation")
        self.assertEqual(self.run_script(), "0.20.7")

    def test_initial_push_only_considers_the_pushed_commit(self):
        self.commit("major: an old breaking change")
        self.commit("moderate: add a feature")
        self.before = "0" * 40
        self.assertEqual(self.run_script(), "0.21.0")

    def test_unknown_previous_revision_falls_back_to_the_pushed_commit(self):
        self.commit("minor: tidy up")
        self.before = "b" * 40
        self.assertEqual(self.run_script(), "0.20.8")


if __name__ == "__main__":
    unittest.main()
