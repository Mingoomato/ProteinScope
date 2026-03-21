"""ProteinScope — Comprehensive Protein Information Aggregator.

CLI entry point. Run:
    python main.py --help
    python main.py P68871 --format markdown --output ./reports/hemoglobin
"""

from __future__ import annotations

import sys
from pathlib import Path

import click
from dotenv import load_dotenv

# Load .env from project root (one level up from this file if run from within proteinscope/)
load_dotenv(Path(__file__).parent / ".env", override=True)

# Ensure the package root is on the Python path
_ROOT = str(Path(__file__).parent)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from core.query_engine import run_query


@click.command()
@click.argument("query")  # UniProt accession, protein name, or gene symbol
@click.option("--organism", default=None, help="Organism filter, e.g. 'Homo sapiens'")
@click.option(
    "--format", "fmt",
    default="markdown",
    type=click.Choice(["json", "markdown", "pdf"], case_sensitive=False),
    help="Output format (default: markdown)",
    show_default=True,
)
@click.option(
    "--output", "-o",
    default="./report",
    help="Output file path without extension (default: ./report)",
    show_default=True,
)
@click.option(
    "--cross-species",
    is_flag=True,
    default=False,
    help="Include cross-species interaction domain comparison (slower)",
)
@click.option(
    "--no-cache",
    is_flag=True,
    default=False,
    help="Bypass local cache and always fetch fresh data",
)
@click.option(
    "--threshold",
    default=80.0,
    type=float,
    help="Identity % threshold for cross-species compatibility (default: 80)",
    show_default=True,
)
@click.option(
    "--no-pathways",
    is_flag=True,
    default=False,
    help="Skip all pathway fetching (faster, sequence/structure only)",
)
@click.option(
    "--pathways-only",
    is_flag=True,
    default=False,
    help="Fetch only pathway data, skip sequence/structure sections",
)
@click.option(
    "--pgx-only",
    is_flag=True,
    default=False,
    help="Fetch only PharmGKB drug/PGx data",
)
@click.option(
    "--min-pgx-evidence",
    default="4",
    show_default=True,
    help="Minimum PharmGKB evidence level (1A/1B/2A/2B/3/4). Default 4 = include all.",
)
@click.option(
    "--interaction-score",
    default=0.4,
    type=float,
    show_default=True,
    help="Minimum STRING confidence score 0.0–1.0 (default: 0.4)",
)
@click.option(
    "--no-images",
    is_flag=True,
    default=False,
    help="Skip pathway diagram image downloads (links only in report)",
)
@click.option(
    "--ai-summary",
    is_flag=True,
    default=False,
    help="Generate AI plain-English summary (section 15) using Gemini 2.5 Pro. "
         "Requires GEMINI_API_KEY in .env.",
)
@click.option(
    "--no-annotate",
    is_flag=True,
    default=False,
    help="Skip diagram annotation (download plain images only, no callouts)",
)
@click.option(
    "--max-callouts",
    default=10,
    type=int,
    show_default=True,
    help="Maximum number of annotated callouts per pathway diagram",
)
@click.option(
    "--no-ai-polish",
    is_flag=True,
    default=False,
    help="Skip Gemini AI polishing of annotation comments (use rule-based text only)",
)
def search(
    query: str,
    organism,
    fmt: str,
    output: str,
    cross_species: bool,
    no_cache: bool,
    threshold: float,
    no_pathways: bool,
    pathways_only: bool,
    pgx_only: bool,
    min_pgx_evidence: str,
    interaction_score: float,
    no_images: bool,
    ai_summary: bool,
    no_annotate: bool,
    max_callouts: int,
    no_ai_polish: bool,
):
    """Search for comprehensive protein information.

    QUERY can be a UniProt accession (e.g. P68871), a gene symbol (e.g. HBB),
    or a protein name (e.g. "hemoglobin beta").

    \b
    Examples:
      python main.py P68871 --format markdown --output ./reports/hemoglobin
      python main.py EGFR --organism "Homo sapiens" --format pdf --output ./reports/egfr
      python main.py "insulin receptor" --format pdf --cross-species
      python main.py P04637 --format json --output ./reports/tp53
      python main.py CYP2C9 --pathways-only --format markdown
      python main.py P13569 --no-pathways --format pdf
      python main.py EGFR --min-pgx-evidence 1A --no-images
    """
    include_pathways = not no_pathways
    try:
        run_query(
            query=query,
            organism=organism,
            fmt=fmt,
            output=output,
            cross_species=cross_species,
            no_cache=no_cache,
            include_pathways=include_pathways,
            no_images=no_images,
            min_pgx_evidence=min_pgx_evidence,
            interaction_score=interaction_score,
            ai_summary=ai_summary,
            no_annotate=no_annotate,
            max_callouts=max_callouts,
            no_ai_polish=no_ai_polish,
        )
    except Exception as exc:
        click.echo(f"Error: {exc}", err=True)
        raise SystemExit(1)


@click.command()
@click.option("--host", default="127.0.0.1", show_default=True, help="Host to bind to")
@click.option("--port", default=8000, show_default=True, type=int, help="Port to listen on")
def serve(host: str, port: int):
    """Launch the ProteinScope web interface.

    \b
    Examples:
      python main.py serve
      python main.py serve --port 8080
      python main.py serve --host 0.0.0.0 --port 8000
    """
    try:
        import uvicorn
    except ImportError:
        click.echo("Error: uvicorn not installed. Run: pip install uvicorn", err=True)
        raise SystemExit(1)
    click.echo(f"Starting ProteinScope web UI at http://{host}:{port}")
    from web.app import app as web_app
    uvicorn.run(web_app, host=host, port=port)


@click.group()
def cli():
    pass


cli.add_command(search)
cli.add_command(serve)


if __name__ == "__main__":
    cli()
