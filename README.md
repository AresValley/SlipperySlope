# SlipperySlope

A simple slippage parameter calcualtor (only LISP implemented so far...)

## How to run

The input coordinate file `[InputXYZ.xyz]` must be in a standard XMol format.

```
python3 SlipperySlope.py [InputXYZ.xyz]
```

- **Index for the ring**: Identifies the idexes of the atoms involved in the reference ring. Must be a ***list of integers*** with a space as separator.
- **Index for the M atom**: This is the index for the M atom. Must be a single ***integer***.