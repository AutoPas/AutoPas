# Build stage with Spack pre-installed and ready to be used
FROM registry.jsc.fz-juelich.de/slms/loadbalancing/spack-alpine as builder


# What we want to install and how we want to install it
# is specified in a manifest file (spack.yaml)
RUN mkdir /opt/spack-environment \
&&  (echo "spack:" \
&&   echo "  specs:" \
&&   echo "  - cmake" \
&&   echo "  - openmpi" \
&&   echo "  - graphviz +expat +libgd" \
&&   echo "  - boost" \
&&   echo "  - doxygen +graphviz" \
&&   echo "  - py-sphinx" \
&&   echo "  - py-sphinx-rtd-theme" \
&&   echo "  - py-breathe" \
&&   echo "  - py-recommonmark" \
&&   echo "  - py-commonmark" \
&&   echo "  - hugo" \
&&   echo "  packages:" \
&&   echo "    all:" \
&&   echo "      target:" \
&&   echo "      - x86_64" \
&&   echo "      compiler:" \
&&   echo "      - gcc@10" \
&&   echo "  concretization: together" \
&&   echo "  config:" \
&&   echo "    install_tree: /opt/software" \
&&   echo "  view: /opt/view") > /opt/spack-environment/spack.yaml

# Install the software, remove unnecessary deps
RUN  cd /opt/spack-environment && spack env activate . && export SPACKMON_USER=$(cat /run/secrets/su) && export SPACKMON_TOKEN=$(cat /run/secrets/st) && spack install --fail-fast && spack gc -y

# Strip all the binaries
RUN find -L /opt/view/* -type f -exec readlink -f '{}' \; | \
    xargs file -i | \
    grep 'charset=binary' | \
    grep 'x-executable\|x-archive\|x-sharedlib' | \
    awk -F: '{print $1}' | xargs strip -s

# Modifications to the environment that are necessary to run
RUN cd /opt/spack-environment && \
    spack env activate --sh -d . >> /etc/profile.d/z10_spack_environment.sh


# Bare OS image to run the installed executables
FROM registry.jsc.fz-juelich.de/slms/loadbalancing/spack-alpine

COPY --from=builder /opt/spack-environment /opt/spack-environment
COPY --from=builder /opt/software /opt/software
COPY --from=builder /opt/view /opt/view
COPY --from=builder /etc/profile.d/z10_spack_environment.sh /etc/profile.d/z10_spack_environment.sh


RUN source /etc/profile && \
wget https://www.vtk.org/files/release/9.0/VTK-9.0.3.tar.gz -O- | tar xz && \
mv VTK-9.0.3 vtk && \
CC=gcc CXX=g++ cmake -S vtk -B vtk_build \
-DCMAKE_BUILD_TYPE=Release \
-DINSTALL_CMAKE_PACKAGE_MODULE=ON \
-DVTK_PYTHON_VERSION=3 \
-DVTK_USE_MPI=ON \
-DVTK_GROUP_ENABLE_MPI=YES \
-DVTK_GROUP_ENABLE_Imaging=DONT_WANT \
-DVTK_GROUP_ENABLE_Qt=DONT_WANT \
-DVTK_GROUP_ENABLE_Rendering=DONT_WANT \
-DVTK_GROUP_ENABLE_StandAlone=WANT \
-DVTK_GROUP_ENABLE_Views=DONT_WANT \
-DVTK_GROUP_ENABLE_Web=DONT_WANT \
-DVTK_ENABLE_LOGGING=OFF && \
cmake --build vtk_build -j16 && \
find -L vtk_build -type f -exec readlink -f '{}' \; | \
xargs file -i | \
grep 'charset=binary' | \
grep 'x-executable\|x-archive\|x-sharedlib' | \
awk -F: '{print $1}' | xargs strip -s && \
cmake --install vtk_build && \
rm -rf vtk vtk_build
RUN sed -i '1s:.*:#!/usr/bin/python3:' /opt/view/bin/sphinx-build /opt/view/bin/sphinx-apidoc /opt/view/bin/sphinx-autogen /opt/view/bin/sphinx-quickstart
RUN source /etc/profile && mkdir /usr/share/fonts && wget https://download.gimp.org/pub/gimp/fonts/urw-fonts.tar.gz -O- | tar xz -C /usr/share/fonts && fc-cache


ENTRYPOINT ["/bin/bash", "--rcfile", "/etc/profile", "-l", "-c"]
