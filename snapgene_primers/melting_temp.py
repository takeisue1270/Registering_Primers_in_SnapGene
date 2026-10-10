"""A deliberately simple melting-temperature estimate.

This value is only ever written into the informational ``meltingTemperature``
attribute that SnapGene displays in its primer list -- it plays no part in
deciding where a primer binds, and is never checked by the verification
step. It is a standard two-regime approximation (Wallace rule for short
oligos, GC%-based formula for longer ones); SnapGene's own internal
nearest-neighbor calculation will differ slightly and that is fine.
"""

from __future__ import annotations


def estimate_tm(sequence: str) -> int:
    seq = sequence.upper()
    length = len(seq)
    if length == 0:
        return 0
    gc = sum(1 for base in seq if base in "GC")
    at = sum(1 for base in seq if base in "AT")
    if length < 14:
        tm = 2 * at + 4 * gc
    else:
        tm = 64.9 + 41 * (gc - 16.4) / length
    return round(tm)
