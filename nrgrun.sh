#!/bin/bash
if [[ -f "data" ]]; then
  if head -n3 "data" | grep -q COMPLEX
  then
    EXEC=nrgcmpl
  else
    EXEC=nrg
  fi
  ${EXEC} "$@" 2>log2 | tee log
else
  echo Input file \"data\" does not exist. Stopping.
fi
