#!/usr/bin/env bash

# https://betterdev.blog/minimal-safe-bash-script-template/
set -Eeuo pipefail

# Set stdout as default so the script can be ran on the commandline.
STEP_SUMMARY="${GITHUB_STEP_SUMMARY:-/dev/stdout}"

echo "# Running Notebooks Tests" >> $STEP_SUMMARY

if [[ "x${1}" == "x--clone-irdb" ]] ; then
  # Cloning IRDB
  if [[ ! -e inst_pkgs ]] ; then
    echo "_Cloning IRDB_" >> $STEP_SUMMARY
    git clone https://github.com/AstarVienna/irdb.git inst_pkgs
  fi

  echo "## Patch notebooks" >> $STEP_SUMMARY
  # https://github.com/koalaman/shellcheck/wiki/SC2044
  find ./docs -iname "*.ipynb" -printf '%h\0' | sort -z | uniq -z | while IFS= read -r -d '' dirnotebooks; do
    echo "${dirnotebooks}"
    echo "- ${dirnotebooks}" >> $STEP_SUMMARY
    pushd "${dirnotebooks}" || exit 1
      # Comment out any download_package[s] in the notebooks.
      sed -i -E 's|"(.*\.download_package)|"#\1|g' -- *.ipynb
      # Comment out explicitly setting the local_packages_path to somewhere.
      sed -i -E 's|"(.*__config__\[\\"!SIM.file.local_packages_path)|"#\1|g' -- *.ipynb
    popd || exit 1
  done
fi


echo "## Notebooks tested" >> $STEP_SUMMARY
# https://github.com/koalaman/shellcheck/wiki/SC2044
find ./docs -iname "*.ipynb" -print0 | while IFS= read -r -d '' fnnotebook
do
  echo "Testing ${fnnotebook} ..."
  fnpy="${fnnotebook%.ipynb}.py"

  # Convert .ipynb file to .py.
  poetry run jupytext --to py "${fnnotebook}"

  # Run the python script and quit on first error.
  poetry run python "${fnpy}" || exit 1
  echo "- ${fnnotebook}" >> $STEP_SUMMARY

  # Delete generated files if --delete is specified.
  # By default do not delete any files.
  # The delete functionality is intended to make it easy for developers to test
  # all the notebooks on their own machine.
  if [ "x$1" = "x--delete" ]
  then
    echo "Removing ${fnpy}"
    rm "${fnpy}"
  fi

done
