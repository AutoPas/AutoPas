# make sure profile is really loaded
# source /etc/profile

# create a badge in current directory
function create_badge {
    if [ -x "$(command -v anybadge)" ]; then
       rm -f "$1"
       [[ ${DETAILS} ]] && echo "BADGE: creating badge with label=$2, value=$3, $4"
       anybadge --file="$1" --label="$2" --value="$3" $4
       rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi
       stat --printf="SIZE \"$1\" -> %s byte\n" $1
    else
       echo "WARNING: Cannot find command anybadge."
    fi
}

# check if command exists & perhaps more
function check_command {
    if ! [ -x "$(command -v $1)" ]; then
       echo "ERROR: Cannot find command $1 ." >&2
       exit 1
    fi
}

# get max. threads available
function set_THREADS {
    if [ -z "$THREADS" ]; then
       [[ ${DETAILS} ]] && echo "The 'THREADS' environment variable is not set, default value = 1."
       THREADS=1
    fi
}

# load MPI
function load_MPI {
#    source /etc/profile.d/modules.sh
#    module add mpi/openmpi-x86_64
    echo "load MPI libs"
}

function find_ALL {
	ALL_ROOT_DIR="$CI_SCRIPT_PATH/.."
}

function find_VTK {
	#VTK_DIR="$CI_SCRIPT_PATH/../../vtk_bin" #local machine
	VTK_DIR=/usr/local
}

# push badge to tmp-branch and exit
function pushbadge_exit {
    local CP_FILE=$1
    local EXIT_CODE=$2
  
    echo " cp -af ${CP_FILE} ${BASE}/badges/badge_${CI_JOB_NAME}_${CI_COMMIT_REF_NAME}.svg"

#    if [[ ${BADGE_TOKEN} ]]; then
      if [[  ${CI_REPOSITORY_URL} ]] && [[ ${CI_COMMIT_REF_NAME} ]] && [[ ${CI_JOB_NAME} ]] ; then

        [[ ${DETAILS} ]] && echo "CI_REPOSITORY_URL  = ${CI_REPOSITORY_URL}"
        [[ ${DETAILS} ]] && echo "CI_COMMIT_REF_NAME = ${CI_COMMIT_REF_NAME}"
        [[ ${DETAILS} ]] && echo "CI_JOB_NAME        = ${CI_JOB_NAME}"
        [[ ${DETAILS} ]] && echo "File to add        = ${CP_FILE}"
        [[ ${DETAILS} ]] && echo "Exit code          = ${EXIT_CODE}"

        cp -af ${CP_FILE} ${BASE}/badges/badge_${CI_JOB_NAME}_${CI_COMMIT_REF_NAME}.svg
        
#        ${CI_SCRIPT_PATH}/git_fileadd.sh \
#          ${CI_REPOSITORY_URL} \
#          ${BADGE_TOKEN} \
#          tmp_${CI_COMMIT_REF_NAME} \
#          ${CP_FILE} \
#          badge_${CI_JOB_NAME}_${CI_COMMIT_REF_NAME}.svg
      fi
#    fi

    exit ${EXIT_CODE}
}

# read command line arguments
POSITIONAL=()
while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
    -t|--token)
    BADGE_TOKEN="$2"
    shift # past argument
    shift # past value
    ;;
    -f|--badge-filename)
    BADGE_FILENAME="$2"
    shift # past argument
    shift # past value
    ;;
    --badge-only)
    BADGE_ONLY=YES
    shift # past argument
    ;;
    -s|--suppressions)
    VALGRIND_SUPPRESSIONS="$2"
    shift # past argument
    shift # past value
    ;;
    --details)
    DETAILS=YES
    shift # past argument
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
  esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

# set -x # print all executed commands to the terminal

# set var THREADS depending on no. processors
set_THREADS

# set basic vars
BASE=$(pwd)
CMAKE=cmake
CI_SCRIPT_PATH=$(cd "$(dirname "$0")"; pwd -P)
if [[ -z ${BADGE_FILENAME} ]]; then
  BADGE_FILENAME="${BASE}/badge.svg"
fi
if [ "${VALGRIND_SUPPRESSIONS}" == "${VALGRIND_SUPPRESSIONS#/}" ]; then
  VALGRIND_SUPPRESSIONS="${BASE}/${VALGRIND_SUPPRESSIONS}"
fi

