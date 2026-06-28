# Hybrid SQLite and Zarr Store

Each Store Release uses a small root manifest, an embedded SQLite database for relational metadata and lookup structures, and a Zarr hierarchy for numerical association arrays. Both Dense and Ragged layouts use this envelope: this preserves SQLite's strengths for metadata and indexing while giving all association payloads the same compressed, chunked numerical substrate instead of forcing either relational data into Zarr or large matrices into SQLite.
