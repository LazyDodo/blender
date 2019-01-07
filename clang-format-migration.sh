#!/bin/bash

FILES=$(
	git ls-tree -r HEAD \
		intern/clog/ \
		intern/ghost/ \
		intern/guardedalloc/ \
		intern/string/ \
		source/ \
		--name-only |
		egrep \\.\(c\|cc\|cpp\|cxx\|h\|hh\|hpp\|hxx\|m\|mm\|osl\|glsl\)$
	 )

# First expand tabs
SCRIPT=$(cat <<EOF
import sys
for f in sys.stdin.read().splitlines():
	print(f)
	with open(f, 'r', encoding="utf-8") as fh:
		data = fh.read()
	data = data.expandtabs(4)
	with open(f, 'w', encoding="utf-8") as fh:
		fh.write(data)

EOF
)
python -c "$SCRIPT" <<< $FILES

# Run clang-format.

# Limit to 2gib (some files can use many gig's of memory).
ulimit -Sv 2000000
xargs clang-format -verbose -i <<< $FILES
