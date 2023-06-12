#!/usr/bin/env bash

if [[ "${1}x" == "--clone-irdbx" ]] ; then
  # Cloning IRDB
  if [[ ! -e irdb ]] ; then
    git clone https://github.com/AstarVienna/irdb.git
  fi

  # https://github.com/koalaman/shellcheck/wiki/SC2044
  find . -iname "*.ipynb" -printf '%h\0' | sort -z | uniq -z | while IFS= read -r -d '' dirnotebooks; do
    echo "${dirnotebooks}"
    dirinstpkgs="${dirnotebooks}/inst_pkgs"
    if [[ (! -e ./docs/source/examples/inst_pkgs) && (! -L ./docs/source/examples/inst_pkgs) ]] ; then
      echo "Cretaing symlink to irdb: ${dirinstpkgs}"
      ln -s irdb "${dirinstpkgs}"
    else
      echo "Dericetory exists, not creating symlink: ${dirinstpkgs}"
    fi
  done
fi


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
