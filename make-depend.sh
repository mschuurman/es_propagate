#!/bin/bash
#
#  Generate a simple Fortran-90 dependency file for module includes
#
for file in "$@" ; do
  base="$(echo "$file" | sed -e 's/.f90$//')"
  echo -n "${base}.o: ${file} "
  awk '/^[ \t]*module[ \t]+/{self[$2]=$2}
       /^[ \t]*use[ \t]+/{if (!($2 in self)){printf "%s.o\n", $2}}' "${file}" | sort -u | tr -s '\12' ' '
  echo "" ; echo ""
done
