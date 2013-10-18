#!/usr/bin/env python
# trapfep.py - integrate compute fep results using the trapezoidal rule

import sys, math

if len(sys.argv) < 3:
    print "usage: trapfep.py temperature hderiv"
    print "output of compute fep read from stdin"
    sys.exit()

rt = 0.008314 * float(sys.argv[1])
hderiv = float(sys.argv[2])

line = sys.stdin.readline()
while line.startswith("#"):
    line = sys.stdin.readline()
tok = line.split()
if len(tok) == 4:
    v = float(tok[3])
else:
    v = 1
lo = -rt * math.log(float(tok[2]) / v)
    
i = 1
sum = 0.0
for line in sys.stdin:
    tok = line.split()
    if len(tok) == 4:
        v = float(tok[3])
    else:
        v = 1
    hi = - rt * math.log(float(tok[2]) / v)
    sum += (hi + lo) / (2 * hderiv)
    lo = hi
    i += 1

print sum / i
