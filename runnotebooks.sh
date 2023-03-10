#!/usr/bin/env bash

# https://github.com/koalaman/shellcheck/wiki/SC2044
find . -iname "*.ipynb" -print0 | while IFS= read -r -d '' fnnotebook
do
  echo "Testing ${fnnotebook}"
  fnpy="${fnnotebook%.ipynb}.py"

  # Convert .ipynb file to .py.
  jupytext --to py "${fnnotebook}"

  # Run the python script and quit on first error.
  python "${fnpy}" || exit 1

  # Delete generated files if --delete is specified.
  # By default do not delete any files.
  if [ "x$1" = "x--delete" ]
  then
    rm "${fnpy}"
  fi

done
