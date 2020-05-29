History
=======

0.1.0 (2019-04-10)
---------------------

-   First release on PyPI.

1.0.0 (2020-04-28)
---------------------

-   Sencha release on PyPI.


1.1.0dev
--------

- Improved indexing to prevent user errors of creating too small of khmer Nodegraph (aka "bloom filter") indices for protein-coding translations
- Added check for "occupancy" of the observed k-mers in the reference proteome as a fraction of the total theoretical k-mer, with the desired fraction being 1e-4 or less. Raises an error and warns the user that the k-mer size is too small.
