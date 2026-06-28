# Store-local variant indices

Each Store Release contains its own union variant table and assigns compact local indices to its rows. Indices need not remain stable across Stores or releases: they exist to connect association arrays to that release's variant table, avoiding the external coordination and dependency of a global numeric variant dictionary.

Reference-Completed releases record Reference Panel Membership as Store variant metadata so query results can distinguish LD-panel variants from observed off-panel variants without adding per-association fields.
