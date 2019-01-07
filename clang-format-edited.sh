#!/bin/bash

FILES=$(
	git ls-tree -r HEAD \
		intern/guardedalloc/ \
		intern/string/ \
		source/ \
		--name-only |
		egrep \\.\(c\|cc\|cpp\|cxx\|h\|hh\|hpp\|hxx\|m\|mm\|osl\|glsl\)$
	 )

if [ -z "$FILES" ]; then
	echo "Nothing to clang-format, exiting!"
	exit 0
fi

# Run clang-format.

# Limit to 2gib (some files can use many gig's of memory).
ulimit -Sv 2000000
xargs clang-format -verbose -i <<< $FILES
