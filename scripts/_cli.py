"""Shared helpers for repository CLI entrypoints."""

from __future__ import annotations

import argparse
from collections.abc import Sequence

Argv = Sequence[str] | None


class CliHelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    """Only show meaningful defaults in `--help` output."""

    def _get_help_string(self, action: argparse.Action) -> str:
        help_text = action.help or ""
        if "%(default)" in help_text:
            return help_text
        if action.default in (None, argparse.SUPPRESS) or action.required:
            return help_text
        return f"{help_text} (default: %(default)s)"


def new_parser(description: str) -> argparse.ArgumentParser:
    return argparse.ArgumentParser(
        allow_abbrev=False,
        description=description,
        formatter_class=CliHelpFormatter,
    )


__all__ = ["Argv", "new_parser"]
