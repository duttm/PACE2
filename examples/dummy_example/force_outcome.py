#!/usr/bin/env python

import sys
x=sys.argv[1]

f = open("out.txt", "a")
f.write(str(x))
f.close()
