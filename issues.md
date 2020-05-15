# Issues with the stochastic model

## Pattern Beginnings

The effects of trans-elements binding to the primary target are local.
For this reaosn, it makes sense to have the `effect` be in terms of a position __relative__ to the pattern.

It's unclear to me whether patterns always unambiguously have a position where they match.

## Pattern interation
Currently we are iterating through each match of a pattern.
I worry this will create too many matches.
For instance consider the G-run, which matches any sequence of at least three Gs.

The sequence "GGGAGGGGGGGGGG" has 37 different matches of G-runs.

It seems likely that not all of these are relevant.

In contrast, Armando's code associated every position in the sequence with at most one match.
This has its own problems:
-  Probably as an oversight, the code is such that if multiple patterns match in a position, only the last one that matches will make a difference.
- Without further constraints on what a match of a pattern is, it is easy to see how this could lead to unintutive results
- Overlapping patterns can occur, but only in a manner that is dependent on the left-to-right order.  Actually, I'm not sure if this is true

## What binds next?

In terms of the binding process there are two questions:
- which translation will bind next
- where will it bind?

## Unbinding

- Do we need to model unbinding?
- After a split, should everything be unbound?

