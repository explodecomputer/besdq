# One reference assembly per Store Release

Every Store Release declares exactly one reference assembly. ALID remains the compact `chr:pos:A1:A2` identity within a Store, while federation and cross-Store matching qualify it with the assembly; liftover is an ingestion concern recorded in provenance rather than permitting mixed coordinate systems inside a Store.
