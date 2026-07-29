#!/usr/bin/env python3
"""Compare separation across modes and write best_mode_report.json."""

from __future__ import annotations

import json
from pathlib import Path

import click
from rich.console import Console

console = Console()


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--today-base",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.option("--target-ez-cutoff", default=3.0, show_default=True, type=float)
def main(today_base: Path, target_ez_cutoff: float) -> None:
    modes = []
    candidates = [
        ("fixed_combo", today_base / "fixed_combo" / "ref_free_ezscore" / "aggregate_summary.json"),
        (
            "filtered_combos",
            today_base / "filtered_combos" / "ref_free_ezscore" / "aggregate_summary.json",
        ),
        (
            "filtered_best_subset",
            today_base
            / "filtered_combos"
            / "ref_free_ezscore"
            / "ez_pair_subset_search.json",
        ),
    ]
    for name, path in candidates:
        if not path.is_file():
            console.print(f"[yellow]skip[/yellow] {name}: missing {path}")
            continue
        data = json.loads(path.read_text())
        if name == "filtered_best_subset":
            sep = data.get("separation_eval_best", {}).get("sep")
            modes.append(
                {
                    "mode": name,
                    "sep_eval_ez": sep,
                    "n_pairs": data.get("best", {}).get("n_pairs"),
                    "source": str(path),
                    "detail": data.get("best", {}),
                }
            )
        else:
            ez = data.get("separation_eval", {}).get("ezscore", {})
            # keys may be str after json
            key = target_ez_cutoff
            if str(key) in ez:
                sep = ez[str(key)].get("sep")
            elif key in ez:
                sep = ez[key].get("sep")
            else:
                # try float keys serialized
                sep = None
                for k, v in ez.items():
                    if abs(float(k) - target_ez_cutoff) < 1e-9:
                        sep = v.get("sep")
                        break
            modes.append(
                {
                    "mode": name,
                    "sep_eval_ez": sep,
                    "n_pairs": data.get("n_ez_combos"),
                    "n_repeats": data.get("total_repeats"),
                    "source": str(path),
                }
            )

    ranked = sorted(
        modes,
        key=lambda m: (-(m["sep_eval_ez"] if m["sep_eval_ez"] is not None else -1.0)),
    )
    report = {
        "target_ez_cutoff": target_ez_cutoff,
        "modes": ranked,
        "best_mode": ranked[0] if ranked else None,
    }
    out = today_base / "best_mode_report.json"
    out.write_text(json.dumps(report, indent=2) + "\n")
    console.print(f"[green]OK[/green] Wrote {out}")
    if ranked:
        b = ranked[0]
        console.print(
            f"[bold]Best mode[/bold]: {b['mode']} "
            f"sep_eval_ez={b['sep_eval_ez']:.4f}"
            if b["sep_eval_ez"] is not None
            else f"[bold]Best mode[/bold]: {b['mode']}"
        )
        for m in ranked:
            console.print(
                f"  {m['mode']}: sep={m.get('sep_eval_ez')} n_pairs={m.get('n_pairs')}"
            )


if __name__ == "__main__":
    main()
