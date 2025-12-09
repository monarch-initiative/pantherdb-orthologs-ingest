"""Command line interface for pantherdb-orthologs-ingest."""
import logging
from pathlib import Path

import typer
from kghub_downloader.download_utils import download_from_yaml
from kghub_downloader.model import DownloadOptions
from koza.runner import KozaRunner
from koza.model.formats import OutputFormat

app = typer.Typer()
logger = logging.getLogger(__name__)


@app.callback()
def callback(
    version: bool = typer.Option(False, "--version", is_eager=True),
):
    """pantherdb-orthologs-ingest CLI."""
    if version:
        from pantherdb_orthologs_ingest import __version__

        typer.echo(f"pantherdb-orthologs-ingest version: {__version__}")
        raise typer.Exit()


@app.command()
def download(force: bool = typer.Option(False, help="Force download of data, even if it exists")):
    """Download data for pantherdb-orthologs-ingest."""
    typer.echo("Downloading data for pantherdb-orthologs-ingest...")
    download_config = Path(__file__).parent / "download.yaml"
    download_options = DownloadOptions(ignore_cache=force)
    download_from_yaml(yaml_file=download_config, output_dir=".", download_options=download_options)


@app.command()
def transform(
    output_dir: str = typer.Option("output", help="Output directory for transformed data"),
    row_limit: int = typer.Option(0, help="Number of rows to process"),
    verbose: bool = typer.Option(False, help="Whether to be verbose"),
):
    """Run the Koza transform for pantherdb-orthologs-ingest."""
    typer.echo("Transforming data for pantherdb-orthologs-ingest...")
    transform_config = Path(__file__).parent / "transform.yaml"
    config, runner = KozaRunner.from_config_file(
        str(transform_config),
        output_dir=output_dir,
        output_format=OutputFormat.tsv,
        row_limit=row_limit,
    )
    runner.run()
    

if __name__ == "__main__":
    app()
