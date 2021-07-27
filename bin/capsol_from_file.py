#!/usr/bin/env python
import capsol
import sys

if len(sys.argv) == 1:
    fname = 'capsol.in'
else:
    fname = sys.argv[1]

print("Capsol Py")
print(f"Running Emily's new capsol function on {fname}")
print("Exiting...")