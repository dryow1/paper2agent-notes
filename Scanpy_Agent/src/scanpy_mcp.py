"""Scanpy MCP server: tools from the "Preprocessing and clustering" tutorial.

Run with: python scanpy_mcp.py  (stdio transport)
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from fastmcp import FastMCP  # noqa: E402

from tools.clustering import clustering_mcp  # noqa: E402

mcp = FastMCP(name="scanpy")
mcp.mount(clustering_mcp)

if __name__ == "__main__":
    mcp.run()
