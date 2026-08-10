"""
HydroLib command-line interface.

Provides CLI commands for validation and benchmarking of the
Bulletin 17C implementation.
"""

from __future__ import annotations

import click


@click.group()
def cli() -> None:
    """HydroLib - Hydrologic frequency analysis tools."""
    pass


@cli.command()
def validate() -> None:
    """Run validation benchmarks against expected PeakfqSA results."""
    from hydrolib.validation.benchmarks import print_benchmark_report, run_all_benchmarks

    click.echo("Running validation benchmarks...")
    results = run_all_benchmarks()
    print_benchmark_report(results)

    n_pass = sum(1 for r in results.values() if r.passed)
    n_total = len(results)
    if n_pass < n_total:
        raise SystemExit(1)


@cli.command()
@click.option("--format", "fmt", type=click.Choice(["text", "json"]), default="text")
def benchmark(fmt: str) -> None:
    """Run benchmarks and generate a report.

    Parameters
    ----------
    fmt : str
        Output format: 'text' or 'json'.
    """
    from hydrolib.validation.benchmarks import run_all_benchmarks
    from hydrolib.validation.reports import generate_json_report, generate_text_report

    click.echo("Running benchmarks...")
    results = run_all_benchmarks()

    if fmt == "json":
        click.echo(generate_json_report(results))
    else:
        click.echo(generate_text_report(results))


@cli.command()
@click.argument("state")
def skew(state: str) -> None:
    """Look up the USGS regional (generalized) skew estimate for STATE.

    STATE may be a full name (e.g. "Vermont") or USPS abbreviation (e.g. "VT").
    States without a dedicated study print the nationwide Bulletin 17B
    fallback instead.
    """
    from hydrolib.regional_skew import get_regional_skew

    try:
        estimate = get_regional_skew(state)
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc

    click.echo(f"State: {state}")
    click.echo(f"Regional skew: {estimate.value:.3f}")
    click.echo(f"Standard error: {estimate.se:.3f}")
    click.echo(f"MSE: {estimate.mse:.3f}")
    click.echo(f"Source: {estimate.source}")
    click.echo(f"Source URL: {estimate.source_url}")
    if estimate.note:
        click.echo(f"Note: {estimate.note}")


if __name__ == "__main__":
    cli()
