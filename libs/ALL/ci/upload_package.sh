#!/bin/bash

PRIVATE_TOKEN=j5NU4Guz_op67gz1SJ-P

# from https://docs.gitlab.com/ee/user/packages/generic_packages/

set -x

ID=SLMS%2Floadbalancing
ID=1403
PKGNAME=ALL-release
PKGVERSION=0.9.0
FNAME=$1

# In CI
# 'curl --header "JOB-TOKEN: $CI_JOB_TOKEN" --upload-file file.txt ${CI_API_V4_URL}/projects/${CI_PROJECT_ID}/packages/generic/my_package/0.0.1/file.txt'

curl --header "PRIVATE-TOKEN: $PRIVATE_TOKEN"\
	--upload-file $(pwd)/$FNAME\
	https://gitlab.version.fz-juelich.de/api/v4/projects/$ID/packages/generic/$PKGNAME/$PKGVERSION/$FNAME
	
