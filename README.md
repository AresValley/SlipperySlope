# SlipperySlope

A simple slippage parameter calcualtor (LISP `[1]` and Basolo Delta)

## How to run

The input coordinate file `[InputXYZ.xyz]` must be in a standard XMol format.

```
python3 SlipperySlope.py [InputXYZ.xyz]
```

- **Index for the ring**: Identifies the idexes of the atoms involved in the reference ring. Must be a ***list of integers*** with a space as separator.
- **Index for the M atom**: This is the index for the M atom. Must be a single ***integer***.

## References

1. [M. Dalla Tiezza, F. M. Bickelhaupt, L. Orian, ***ChemistryOpen*** **2019**, 8, 143](https://onlinelibrary.wiley.com/doi/full/10.1002/open.201800191)