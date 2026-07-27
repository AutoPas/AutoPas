#!/bin/bash

set -x

# Run without any arguments to just create a tarball, which is named
# `all-$GITCOM.tar.bz2` where $GITCOM is either the current tag or commit, if
# it is not tagged.
#
# Run with the argument `upload` this script will also publish the tarball as a
# generic package assuming a CI environment.
#
# Run with the argument `upload-local` it will publish the tarball as a generic
# package assuming it is run from a user's environment. This means the
# environment variable PRIVATE_TOKEN must be set to a private token of the user
# with sufficient access rights to the repository.

GITCOM=$(git describe --exact-match --tags)
GITTAG=YES
if [[ -z $GITCOM ]]
then
	GITCOM=$(git rev-parse --short --verify HEAD)
	GITTAG=NO
	if [[ $? -ne 0 ]]
	then
		echo "GIT ERROR... ($GITCOM)"
		echo "Need git to work"
		exit 1
	fi
fi

echo "Creating tarball for release: $GITCOM"

TARNAME=all-$GITCOM

if [[ $2 != reuse ]]
then

git archive --prefix=$TARNAME/ -o $TARNAME.tar $GITCOM

mkdir ${TARNAME}_build
cmake -S . -B ${TARNAME}_build -DCM_ALL_AUTO_DOC=ON -DCMAKE_INSTALL_PREFIX=$(pwd)/${TARNAME}_install
cmake --build ${TARNAME}_build
cmake --install ${TARNAME}_build

MAINPATH=$(pwd)

mkdir -p $TARNAME/Documentation
cd $TARNAME/Documentation
cp -r "$MAINPATH/${TARNAME}_install/share/doc/ALL/html/" .
mv html Doxygen
cp -r "$MAINPATH/${TARNAME}_install/share/doc/ALL/sphinx/" .
mv sphinx Sphinx
cd ../..

tar cf tmp.tar $TARNAME
tar --concatenate --file=$TARNAME.tar tmp.tar
rm -f tmp.tar

rm -rf ${TARNAME}_build ${TARNAME}_install ${TARNAME}

bzip2 -f $TARNAME.tar
fi

ID=1403
PKGNAME=ALL-release
# $GITCOM is either v1.2.3-something2 or git short ref, i.e. 5 hex chars
PKGVERSION=${GITCOM:1:5}
PUBLISHVERSION=$PKGVERSION
FNAME=$TARNAME.tar.bz2

if [[ ${PKGVERSION:1:1} != . ]]
then
	# Sentinel for untagged and we only have commit hash
	PUBLISHVERSION=0.0.0
fi

if [[ $1 == upload ]]
then
	curl --header "JOB-TOKEN: $CI_JOB_TOKEN"\
		--upload-file $(pwd)/$FNAME\
		${CI_API_V4_URL}/projects/${CI_PROJECT_ID}/packages/generic/$PKGNAME/$PUBLISHVERSION/$FNAME
elif [[ $1 == upload-local ]]
then
	curl --header "PRIVATE-TOKEN: $PRIVATE_TOKEN"\
		--upload-file $(pwd)/$FNAME\
		https://gitlab.jsc.fz-juelich.de/api/v4/projects/$ID/packages/generic/$PKGNAME/$PUBLISHVERSION/$FNAME
fi










