FROM registry.jsc.fz-juelich.de/slms/loadbalancing/spack_v2:latest
LABEL maintainer="Rene Halver <r.halver@fz-juelich.de>"

RUN update-alternatives --force --install /usr/bin/gcc gcc /usr/bin/gcc-10 40 && \
    update-alternatives --force --install /usr/bin/gcc-ar gcc-ar /usr/bin/gcc-ar-10 40 && \
    update-alternatives --force --install /usr/bin/gcc-nm gcc-nm /usr/bin/gcc-nm-10 40 && \
    update-alternatives --force --install /usr/bin/gcc-ranlib gcc-ranlib /usr/bin/gcc-ranlib-10 40 && \
    update-alternatives --force --install /usr/bin/g++ g++ /usr/bin/g++-10 40 && \
    update-alternatives --force --install /usr/bin/gfortran gfortran /usr/bin/gfortran-10 40 && \
    update-alternatives --force --install /usr/bin/gcov gcov /usr/bin/gcov-10 40 && \
    update-alternatives --force --install /usr/bin/gcov-dump gcov-dump /usr/bin/gcov-dump-10 40 && \
    update-alternatives --force --install /usr/bin/gcov-tool gcov-tool /usr/bin/gcov-tool-10 40

ADD packages.yaml /root/.spack/packages.yaml

RUN spack install cmake && \
    /root/spack-rm-confdir3.sh && \
    spack clean -a && \
    spack view --dependencies true symlink --ignore-conflicts /usr/local cmake 

RUN spack install openmpi && \
    /root/spack-rm-confdir3.sh && \
    spack clean -a && \
    spack view --dependencies true symlink /usr/local openmpi

RUN wget https://www.vtk.org/files/release/9.0/VTK-9.0.3.tar.gz -O- | tar xz && \
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

#RUN spack install vtk -osmesa -opengl2 -qt ^openmpi && \
#    /root/spack-rm-confdir3.sh && \
#    spack clean -a && \
#    spack view --dependencies false symlink --ignore-conflicts /usr/local vtk

#RUN spack dependencies graphviz~ghostscript && \
#    spack graph mesa~llvm  && \
#    spack graph graphviz~ghostscript ^osmesa~llvm && \
#    spack install graphviz~ghostscript ^osmesa~llvm && \
#    spack view --dependencies no symlink --ignore-conflicts /usr/local graphviz

RUN apt-get update && \
    apt-get -y install graphviz && \
    apt-get -y install libboost-all-dev

RUN spack install doxygen && \
    spack install py-sphinx && \
    spack install py-breathe && \
    spack install py-recommonmark && \
    spack install py-commonmark && \
    spack install py-pip && \
    spack install hugo && \
    /root/spack-rm-confdir3.sh && \
    spack clean -a && \
    spack view --dependencies no symlink --ignore-conflicts /usr/local doxygen && \
    spack view --dependencies true symlink --ignore-conflicts /usr/local py-sphinx && \
    spack view --dependencies true symlink --ignore-conflicts /usr/local py-pip && \
    spack view --dependencies no symlink --ignore-conflicts /usr/local py-breathe && \
    spack view --dependencies no symlink --ignore-conflicts /usr/local py-recommonmark && \
    spack view --dependencies no symlink --ignore-conflicts /usr/local py-commonmark && \
    spack view --dependencies no symlink --ignore-conflicts /usr/local hugo && \
    pip install sphinx_rtd_theme

