# Retain full normalized long alleles

Canonical variants are trimmed and left-aligned before identity is assigned. Long alleles may be represented compactly by a deterministic hashed ALID and stored through a compressed overflow representation, but the Store retains the complete normalized alleles once per variant so exact export, validation, and harmonisation do not depend on an irreversible hash.
