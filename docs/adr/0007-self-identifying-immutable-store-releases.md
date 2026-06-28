# Self-identifying immutable Store Releases

Published BESDQ Stores are immutable releases whose intrinsic manifest records the stable Store identity, exact release identity, format version, creation time, and provenance. A separate catalogue API owns discovery, deployment locations, and selection of the current release; release identity is also embedded in the Store so downloaded and mirrored copies remain standalone, reproducible, and safe to cache.

An Observed-Only release and its Reference-Completed enhancement share the same Store identity but have different release identities. Completion state is a release property, not a reason to mint a separate Store.
