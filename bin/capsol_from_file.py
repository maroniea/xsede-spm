#!/usr/bin/env python
import capsol.capsol as cp
import sys

if len(sys.argv) == 1:
    fname = 'capsol.in'
else:
    fname = sys.argv[1]

print("Capsol Py")
cp.runnewcapsol(fname)
print("Exiting...")

