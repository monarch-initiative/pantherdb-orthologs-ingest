from pathlib import Path

from kgx.cli.cli_utils import transform as kgx_transform
from loguru import logger


out_report = "output/pantherdb_orthologs.nt.gz"
logger.info("Creating rdf output: {}...".format(out_report))

src_files = []
src_nodes = "output/pantherdb_orthologs_nodes.tsv"
src_edges = "output/pantherdb_orthologs_edges.tsv"

if Path(src_nodes).is_file():
    src_files.append(src_nodes)
if Path(src_edges).is_file():
    src_files.append(src_edges)

kgx_transform(inputs=src_files,
              input_format="tsv",
              stream=True,
              output=out_report,
              output_format="nt",
              output_compression="gz")
