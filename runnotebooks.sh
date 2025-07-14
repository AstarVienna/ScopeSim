#!/usr/bin/env bash

# https://betterdev.blog/minimal-safe-bash-script-template/
set -Eeuo pipefail

# Set stdout as default so the script can be ran on the commandline.
STEP_SUMMARY="${GITHUB_STEP_SUMMARY:-/dev/stdout}"

usage() {
  cat <<EOF
Usage: $(basename "${BASH_SOURCE[0]}") [-h] [-v] [--clone-irdb] [--delete]

Test whether all the notebooks run.

Available options:

-h, --help      Print this help and exit
-v, --verbose   Print script debug info
--clone-irdb    Clone the IRDB instead of downloading it
--delete        Delete generated files
EOF
  exit
}

parse_params() {
  # default values of variables set from params
  DELETE=0
  CLONEIRDB=0

  while :; do
    case "${1-}" in
    -h | --help) usage ;;
    -v | --verbose) set -x ;;
    --no-color) NO_COLOR=1 ;;
    --clone-irdb) CLONEIRDB=1 ;;
    --delete) DELETE=1 ;;
    -?*) die "Unknown option: $1" ;;
    *) break ;;
    esac
    shift
  done

  args=("$@")

  return 0
}

msg() {
  echo >&2 -e "${1-}"
}

die() {
  local msg=$1
  local code=${2-1} # default exit status 1
  msg "$msg"
  exit "$code"
}

parse_params "$@"

echo "# Running Notebooks Tests" >> $STEP_SUMMARY

if [[ "${CLONEIRDB}" == 1 ]] ; then
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
    pushd "${dirnotebooks}"
      # Comment out any download_package[s] in the notebooks.
      sed -i -E 's|"(.*\.download_package)|"#\1|g' -- *.ipynb
      # Comment out explicitly setting the local_packages_path to somewhere.
      sed -i -E 's|"(.*__config__\[\\"!SIM.file.local_packages_path)|"#\1|g' -- *.ipynb
    popd
  done
fi


echo "## Notebooks tested" >> $STEP_SUMMARY
# https://github.com/koalaman/shellcheck/wiki/SC2044
# Sort is necessary because some notebooks are intended to be ran in order.
# E.g. 3_custom_effects.ipynb uses MICADO, which is downloaded by 1_scopesim_intro.ipynb
# Also, sorting makes the test deterministic.
find ./docs -iname "*.ipynb" -print0 | sort -z | while IFS= read -r -d '' fnnotebook
do
  echo "Testing ${fnnotebook} ..."
  fnpy="${fnnotebook%.ipynb}.py"

  # Convert .ipynb file to .py.
  jupytext --to py "${fnnotebook}"

  # Run the python script.
  python "${fnpy}"
  echo "- ${fnnotebook}" >> $STEP_SUMMARY

  # Delete generated files if --delete is specified.
  # By default do not delete any files.
  # The delete functionality is intended to make it easy for developers to test
  # all the notebooks on their own machine.
  if [ "${DELETE}" = 1 ]
  then
    echo "Removing ${fnpy}"
    rm "${fnpy}"
  fi

done
