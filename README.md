# Motif Transition Model

This is a C++ implementation of motif transition model (MTM).

To compile the code
`cd src`
`make`

To run the code
`./MTM inputfile l_max delta output_number`

Format of the inputfile (source node, target node, timestamp)
`u1 v1 t1`
`u2 v2 t2`
`...`

### Reference
```
@inproceedings{10.1145/3534678.3539234,
author={Liu, Penghang and Sarıyüce, Ahmet Erdem},
title = {Using Motif Transitions for Temporal Graph Generation},
year = {2023},
doi = {10.1145/3580305.3599540},
booktitle = {Proceedings of the 29th ACM SIGKDD Conference on Knowledge Discovery and Data Mining}
}

```
