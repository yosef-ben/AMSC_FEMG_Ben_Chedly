# Connectome source files

`budapest_all_20k.graphml` is the public all-subject 20k-fibre graph used to
reconstruct the 83-region topology and median connectivity weights.

`budapest_all_20k_mean_electrical.graphml` was downloaded from the official
Budapest Reference Connectome server with version 3.0, all subjects,
`minOccurrences=5`, zero weight threshold, electrical-connectivity weights,
20k launched fibres, and mean combination mode. The reproducible endpoint is:

```text
https://pitgroup.org/apps/connectome/getgraph.php?format=graphml&version=2&population=0&minOccurrences=5&minStrength=0&combineMode=mean&weightFunction=0&totalFiberNumber=0
```

The mean-mode file is retained for the data audit. After aggregation it is
farther from the connectivity statistics reported by Fornari et al. than the
median source, so it is not silently substituted into benchmark 19.
