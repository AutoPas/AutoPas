FROM gitlab.version.fz-juelich.de:5555/slms/loadbalancing/spack_v2:latest
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

RUN spack install mpich && \
    /root/spack-rm-confdir3.sh && \
    spack clean -a && \
    spack view --dependencies true symlink /usr/local mpich

ADD vtk_shrunk.tar.gz /

RUN mkdir vtk_build && \
    cd vtk_build && \
    CC=gcc CXX=g++ cmake ../vtk \
    $(EXT_LTO_CMFLAGS) \
    -DBUILD_SHARED_LIBS=OFF \
    -DBUILD_TESTING=OFF \
    -DCMAKE_BUILD_TYPE=Release \
    -DVTK_Group_MPI=OFF \
    -DVTK_Group_Rendering=OFF \
    -DVTK_Group_StandAlone=OFF \
    -DVTK_RENDERING_BACKEND=None \
    -DVTK_USE_CXX11_FEATURES=ON \
    -DModule_vtkCommonDataModel=ON \
    -DModule_vtkFiltersProgrammable=ON \
    -DModule_vtkIOParallelXML=ON \
    -DModule_vtkParallelMPI=ON && \
    make -j 4 -C . && \
    make -C . install

#RUN spack install vtk -osmesa -opengl2 -qt ^mpich && \
#    /root/spack-rm-confdir3.sh && \
#    spack clean -a && \
#    spack view --dependencies false symlink --ignore-conflicts /usr/local vtk

#RUN spack dependencies graphviz~ghostscript && \
#    spack graph mesa~llvm  && \
#    spack graph graphviz~ghostscript ^osmesa~llvm && \
#    spack install graphviz~ghostscript ^osmesa~llvm && \
#    spack view --dependencies no symlink --ignore-conflicts /usr/local graphviz

RUN apt-get -y install graphviz && \
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

