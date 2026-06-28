# Local direct queries and remote access through the API

The v1 Python package directly queries downloaded or otherwise local Store Releases; remote querying is provided by the separate web API. Direct HTTP or object-store querying is not part of the v1 Store contract, allowing each release to use an embedded relational index while remaining standalone after download.

Build and query CLI commands operate on explicit Store paths supplied as arguments. Directory naming conventions and multi-store organisation are outside the BESDQ Store contract and belong to a separate catalogue or deployment layer.
