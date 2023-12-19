#!/usr/bin/env bash

echo "# Running Notebooks Tests" >> $GITHUB_STEP_SUMMARY

if [[ "x${1}" == "x--clone-irdb" ]] ; then
  # Cloning IRDB
  if [[ ! -e irdb ]] ; then
    echo "_Cloning IRDB_" >> $GITHUB_STEP_SUMMARY
    git clone https://github.com/AstarVienna/irdb.git
  fi

  echo "## Symlinks" >> $GITHUB_STEP_SUMMARY
  # https://github.com/koalaman/shellcheck/wiki/SC2044
  find . -iname "*.ipynb" -printf '%h\0' | sort -z | uniq -z | while IFS= read -r -d '' dirnotebooks; do
    echo "${dirnotebooks}"
    echo "- ${dirnotebooks}" >> $GITHUB_STEP_SUMMARY
    dirinstpkgs="${dirnotebooks}/inst_pkgs"
    if [[ (! -e ./docs/source/examples/inst_pkgs) && (! -L ./docs/source/examples/inst_pkgs) ]] ; then
      echo "Creating symlink to irdb: ${dirinstpkgs}" >> $GITHUB_STEP_SUMMARY
      ln -s irdb "${dirinstpkgs}"
    else
      echo "Directory exists, not creating symlink: ${dirinstpkgs}" >> $GITHUB_STEP_SUMMARY
    fi

    # Comment out any download_package[s] in the notebooks.
    pushd "${dirnotebooks}" || exit 1
      sed -i -E 's|"(.*\.download_package)|"#\1|g' -- *.ipynb
    popd || exit 1
  done
fi


echo "## Notebooks tested" >> $GITHUB_STEP_SUMMARY
# https://github.com/koalaman/shellcheck/wiki/SC2044
find . -iname "*.ipynb" -print0 | while IFS= read -r -d '' fnnotebook
do
  echo "Testing ${fnnotebook} ..."
  fnpy="${fnnotebook%.ipynb}.py"

  # Convert .ipynb file to .py.
  poetry run jupytext --to py "${fnnotebook}"

  # Run the python script and quit on first error.
  poetry run python "${fnpy}" || exit 1
  echo "- ${fnnotebook}" >> $GITHUB_STEP_SUMMARY

  # Delete generated files if --delete is specified.
  # By default do not delete any files.
  if [ "x$1" = "x--delete" ]
  then
    rm "${fnpy}"
  fi

done
