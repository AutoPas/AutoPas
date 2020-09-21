pipeline {
    agent {
        node {
            //cloud 'kubernetes'
            label 'openshift-autoscale'
        }
    }
    stages {
        stage('setup') {
            steps {
                echo 'Starting AutoPas Pipeline'
            }
        }
        stage("style check") {
            parallel {
                stage("build documentation") {
                    steps {
                        container('autopas-cmake-doxygen-make') {
                            dir("build-doxygen") {
                                sh 'entrypoint.sh ccache -s'
                                sh 'cmake ..'
                                sh 'make doc_doxygen 2>DoxygenWarningLog.txt'
                            }
                        }
                        stash includes: 'build-doxygen/doc_doxygen/html/**', name: 'doxydocs'

                        // get doxygen warnings
                        recordIssues filters: [excludeFile('.*README.*')], tools: [doxygen(id: 'docAutoPas', pattern: 'build-doxygen/DoxygenWarningLog.txt')], failedTotalAll: 1
                    }
                    post {
                        failure {
                            error "warnings in doxygen documentation"
                        }
                    }
                }
                stage("build md-flexible documentation") {
                    steps {
                        container('autopas-cmake-doxygen-make') {
                            dir("build-doxygen-md-flexible") {
                                sh 'entrypoint.sh ccache -s'
                                sh 'cmake ..'
                                sh 'make doc_doxygen_md-flexible 2>DoxygenWarningLog.txt'
                            }
                        }
                        stash includes: 'build-doxygen-md-flexible/doc_doxygen_md-flexible/html/**', name: 'doxydocs_md-flexible'

                        // get doxygen warnings
                        recordIssues filters: [excludeFile('.*README.*')], tools: [doxygen(id: 'docMdFlex', pattern: 'build-doxygen-md-flexible/DoxygenWarningLog.txt')], failedTotalAll: 1
                    }
                    post {
                        failure {
                            error "warnings in doxygen documentation for md-flexible"
                        }
                    }
                }
                stage ("clang and cmake format") {
                    steps {
                        dir("format") {
                            container('autopas-clang6-cmake-ninja-make') {
                                sh "CC=clang CXX=clang++ cmake -G Ninja -DAUTOPAS_OPENMP=ON .."
                                sh "ninja clangformat"
                                sh "ninja cmakeformat"
                            }
                            script {
                                // return 2 if files have been modified by clang-format, 0 otherwise
                                try {
                                    // if files were modified, return 2
                                    sh "git diff --quiet || exit 2"
                                } catch (Exception e) {
                                    // change detected
                                    echo 'clang or cmake format errors detected. please format the code properly. Affected files:'
                                    sh "git status | grep modified"
                                    sh "exit 1"
                                }
                            }
                        }
                    }
                }
                stage ("custom checks") {
                    steps {
                        echo 'Testing src folder'
                        dir("src") {
                            checkCustom()
                        }
                        echo 'Testing tests folder'
                        dir("tests") {
                            checkCustom()
                        }
                        echo 'Testing examples folder'
                        dir("examples") {
                            checkCustom()
                        }
                    }
                }
            }
        }
        stage('build and test') {
            options {
                timeout(time: 6, unit: 'HOURS')
            }
            parallel {
                stage('gpu cloud') {
                    agent { label 'openshift-autoscale-gpu' }
                    steps {
                        container('cuda-10') {
                            dir("build-cuda") {
                                sh "cmake -DAUTOPAS_ENABLE_CUDA=ON -DCCACHE=ON .."
                                sh "entrypoint.sh make -j 4 > buildlog-cuda.txt 2>&1 || (cat buildlog-cuda.txt && exit 1)"
                                sh "./tests/testAutopas/runTests"
                            }
                            dir('build-cuda/examples') {
                                sh "ctest -C checkExamples -j8 --verbose"
                            }
                        }
                    }
                    post {
                        always {
                            recordIssues qualityGates: [[threshold: 1, type: 'TOTAL', unstable: true]], tools: [gcc(id: 'cuda-gcc', pattern: 'build*/buildlog-cuda.txt')]
                        }
                    }
                }
                stage('gpu cloud - clang') {
                    agent { label 'openshift-autoscale-gpu' }
                    steps {
                        container('cuda-10') {
                            dir("build-cuda") {
                                sh "CC=clang CXX=clang++ cmake -DCCACHE=ON -DAUTOPAS_ENABLE_CUDA=ON .."
                                sh "entrypoint.sh make -j 4 > buildlog-cuda-clang.txt 2>&1 || (cat buildlog-cuda-clang.txt && exit 1)"
                                sh "./tests/testAutopas/runTests"
                                // cuda variants of valgrind:
                                sh "cuda-memcheck --tool memcheck ./tests/testAutopas/runTests"
                                sh "cuda-memcheck --tool racecheck ./tests/testAutopas/runTests"
                                sh "cuda-memcheck --tool synccheck ./tests/testAutopas/runTests"
                                // initcheck does not pass, I think because we don't initialize patted values.
                                //sh "cuda-memcheck --tool initcheck ./tests/testAutopas/runTests"
                            }
                            dir('build-cuda/examples') {
                                sh "ctest -C checkExamples -j8 --verbose"
                            }
                        }
                    }
                    post {
                        always {
                            recordIssues qualityGates: [[threshold: 1, type: 'TOTAL', unstable: true]], tools: [clang(id: 'cuda-clang', pattern: 'build*/buildlog-cuda-clang.txt')]
                        }
                    }
                }
                stage("default") {
                    steps {
                        container('autopas-gcc7-cmake-make') {
                            dir("build") {
                                sh "cmake -DCCACHE=ON .."
                                sh "entrypoint.sh make -j 4 > buildlog.txt 2>&1 || (cat buildlog.txt && exit 1)"
                                sh 'env GTEST_OUTPUT="xml:$(pwd)/test.xml" ./tests/testAutopas/runTests'
                            }
                            dir("build/examples") {
                                sh 'ctest -C checkExamples -j8 --verbose'
                            }
                        }
                    }
                }
                stage("gcc openmp") {
                    steps {
                        container('autopas-gcc7-cmake-make') {
                            dir("build-openmp") {
                                sh "cmake -DCCACHE=ON -DAUTOPAS_OPENMP=ON .."
                                sh "entrypoint.sh make -j 4 > buildlog.txt 2>&1 || (cat buildlog.txt && exit 1)"
                                sh './tests/testAutopas/runTests'
                            }
                            dir("build-openmp/examples") {
                                sh 'ctest -C checkExamples -j8 --verbose'
                            }
                        }
                    }
                }
                stage("gcc openmp address-sanitizer") {
                    steps {
                        container('autopas-gcc7-cmake-make') {
                            dir("build-openmp-address-sanitizer") {
                                sh "cmake -DCCACHE=ON -DAUTOPAS_OPENMP=ON -DCMAKE_BUILD_TYPE=Debug -DAUTOPAS_ENABLE_ADDRESS_SANITIZER=ON .."
                                sh "entrypoint.sh make -j 4 > buildlog.txt 2>&1 || (cat buildlog.txt && exit 1)"
                                sh './tests/testAutopas/runTests'
                            }
                        }
                    }
                }
                stage("address sanitizer") {
                    steps {
                        container('autopas-gcc7-cmake-make') {
                            dir("build-addresssanitizer") {
                                sh "cmake -DCCACHE=ON -DCMAKE_BUILD_TYPE=Debug -DAUTOPAS_ENABLE_ADDRESS_SANITIZER=ON .."
                                sh "entrypoint.sh make -j 4 > buildlog.txt 2>&1 || (cat buildlog.txt && exit 1)"
                                sh './tests/testAutopas/runTests'
                            }
                        }
                    }
                }
                stage("address sanitizer release") {
                    steps {
                        container('autopas-gcc7-cmake-make') {
                            dir("build-addresssanitizer-release") {
                                sh "cmake -DCCACHE=ON -DCMAKE_BUILD_TYPE=Release -DAUTOPAS_ENABLE_ADDRESS_SANITIZER=ON .."
                                sh "entrypoint.sh make -j 4 > buildlog.txt 2>&1 || (cat buildlog.txt && exit 1)"
                                sh './tests/testAutopas/runTests'
                            }
                        }
                    }
                }
                stage("clang openmp") {
                    steps {
                        container('autopas-clang6-cmake-ninja-make') {
                            dir("build-clang-ninja-openmp") {
                                sh "CC=clang CXX=clang++ cmake -G Ninja -DCCACHE=ON -DAUTOPAS_OPENMP=ON .."
                                sh "entrypoint.sh ninja -j 4 > buildlog_clang.txt 2>&1 || (cat buildlog_clang.txt && exit 1)"
                                sh './tests/testAutopas/runTests'
                            }
                            dir("build-clang-ninja-openmp/examples") {
                                sh 'ctest -C checkExamples -j8 --verbose'
                            }
                        }
                    }
                }
                stage("clang ninja address sanitizer") {
                    steps {
                        container('autopas-clang6-cmake-ninja-make') {
                            dir("build-clang-ninja-addresssanitizer-debug") {
                                sh "CXXFLAGS=-Wno-pass-failed CC=clang CXX=clang++ cmake -G Ninja -DCCACHE=ON -DCMAKE_BUILD_TYPE=Debug -DAUTOPAS_ENABLE_ADDRESS_SANITIZER=ON .."
                                sh "entrypoint.sh ninja -j 4 > buildlog_clang.txt 2>&1 || (cat buildlog_clang.txt && exit 1)"
                                sh './tests/testAutopas/runTests'
                            }
                        }
                    }
                }
                stage("clang ninja address sanitizer release") {
                    steps {
                        container('autopas-clang6-cmake-ninja-make') {
                            dir("build-clang-ninja-addresssanitizer-release") {
                                sh "CXXFLAGS=-Wno-pass-failed CC=clang CXX=clang++ cmake -G Ninja -DCCACHE=ON -DCMAKE_BUILD_TYPE=Release -DAUTOPAS_ENABLE_ADDRESS_SANITIZER=ON .."
                                sh "entrypoint.sh ninja -j 4 > buildlog_clang.txt 2>&1 || (cat buildlog_clang.txt && exit 1)"
                                sh './tests/testAutopas/runTests'
                            }
                            dir("build-clang-ninja-addresssanitizer-release/examples") {
                                sh 'ctest -C checkExamples -j8 --verbose'
                            }
                        }
                    }
                }
                stage("archer") {
                    steps {
                        container('autopas-llvm-archer') {
                            dir("build-archer") {
                                sh "CXXFLAGS=-Wno-pass-failed CC=clang CXX=clang++ cmake -G Ninja -DAUTOPAS_ENABLE_THREAD_SANITIZER=ON -DAUTOPAS_OPENMP=ON -DCCACHE=ON -DCMAKE_BUILD_TYPE=Release -DAUTOPAS_USE_VECTORIZATION=OFF .."
                                sh "entrypoint.sh ninja -j 4 > buildlog_clang.txt 2>&1 || (cat buildlog_clang.txt && exit 1)"
                                sh 'export TSAN_OPTIONS="ignore_noninstrumented_modules=1" && ctest --verbose -j8'
                            }
                            dir("build-archer/examples") {
                                sh 'export TSAN_OPTIONS="ignore_noninstrumented_modules=1" && ctest -C checkExamples -j8 --verbose'
                            }
                        }
                    }
                }
                stage("intel") {
                    steps {
                        container('autopas-intel18') {
                            dir("build-intel") {
                                sh "bash -i -c 'which icc && CC=`which icc` CXX=`which icpc` cmake -DCCACHE=ON -DAUTOPAS_OPENMP=OFF ..'"
                                sh "bash -i -c 'uid_entrypoint make -j 4 > buildlog_intel.txt 2>&1 || (cat buildlog_intel.txt && exit 1)'"
                                sh "bash -i -c './tests/testAutopas/runTests'"
                            }
                            dir("build-intel/examples") {
                                sh "bash -i -c 'ctest -C checkExamples -j8 --verbose'"
                            }
                        }
                    }
                }
                stage("intel openmp") {
                    steps {
                        container('autopas-intel18') {
                            dir("build-intel-ninja-openmp") {
                                sh "bash -i -c 'which icc && CC=`which icc` CXX=`which icpc` cmake -G Ninja -DCCACHE=ON -DAUTOPAS_OPENMP=ON ..'"
                                sh "bash -i -c 'uid_entrypoint ninja -j 4 > buildlog_intel.txt 2>&1 || (cat buildlog_intel.txt && exit 1)'"
                                sh "bash -i -c './tests/testAutopas/runTests'"
                            }
                            dir("build-intel-ninja-openmp/examples") {
                                sh "bash -i -c 'ctest -C checkExamples -j8 --verbose'"
                            }
                        }
                    }
                }
            }
            post {
                always {
                    recordIssues qualityGates: [[threshold: 1, type: 'TOTAL', unstable: true]], tools: [clang(pattern: 'build*/buildlog_clang.txt'), intel(pattern: 'build*/buildlog_intel.txt'), gcc(pattern: 'build*/buildlog.txt')]
                }
            }
        }
        stage("generate reports") {
            steps {
                // get test results -- mainly to get number of tests
                junit 'build/test.xml'

                // generate coverage
                dir("coverage") {
                    container('autopas-build-code-coverage') {
                        sh "cmake -DAUTOPAS_CODE_COVERAGE=ON -DCCACHE=ON -DCMAKE_BUILD_TYPE=Debug .."
                        sh "entrypoint.sh make AutoPas_cobertura -j 4"
                    }
                    cobertura autoUpdateHealth: false, autoUpdateStability: false, coberturaReportFile: 'coverage.xml', conditionalCoverageTargets: '70, 0, 0', failUnhealthy: false, failUnstable: false, lineCoverageTargets: '80, 0, 0', maxNumberOfBuilds: 0, methodCoverageTargets: '80, 0, 0', onlyStable: false, sourceEncoding: 'ASCII', zoomCoverageChart: false
                }
            }
        }

        stage("publish documentation") {
            when { branch 'master' }
            agent { label 'www_access' }
            steps {
                unstash 'doxydocs'
                dir("build-doxygen") {
                    sh 'touch /import/www/wwwsccs/html/AutoPas/doxygen_doc/master || echo 0'
                    sh 'rm -rf /import/www/wwwsccs/html/AutoPas/doxygen_doc/master || echo 0'
                    sh 'cp -r doc_doxygen/html /import/www/wwwsccs/html/AutoPas/doxygen_doc/master'
                }
            }
        }
        stage("publish md-flexible documentation") {
            when { branch 'master' }
            agent { label 'www_access' }
            steps {
                unstash 'doxydocs_md-flexible'
                dir("build-doxygen-md-flexible") {
                    sh 'touch /import/www/wwwsccs/html/AutoPas/doc_doxygen_md-flexible/master || echo 0'
                    sh 'rm -rf /import/www/wwwsccs/html/AutoPas/doc_doxygen_md-flexible/master || echo 0'
                    sh 'cp -r doc_doxygen_md-flexible/html /import/www/wwwsccs/html/AutoPas/doc_doxygen_md-flexible/master'
                }
            }
        }
    }
    post {
        changed {
            echo "change"
        }
        unstable {
            echo "unstable"
        }
        success {
            echo "success"
        }
        failure {
            echo "failure"
        }
        aborted {
            echo "aborted"
        }
    }
}

void checkCustom() {
    script {
        // check if all header files have a #pragma once
        try {
            // if header files do not contain #pragma once, make build unstable
            sh 'grep -L "#pragma once" -r . | grep -q "\\.h" && exit 2 || exit 0'
        } catch (Exception e) {
            // change detected
            echo 'all header include guards should be implemented using "#pragma once". Affected files:'
            sh 'grep -L "#pragma once" -r . | grep "\\.h"'
            sh "exit 1"
        }

        // check if all files are documented with @file or \file doxygen comments
        try {
            // if .cpp or .h files do not contain a file comment, return 2
            sh "grep '\\\\file\\|\\@file' -Lr . | grep -q '\\.cpp\\|\\.h' && exit 2 || exit 0"
        } catch (Exception e) {
            // change detected
            echo 'all .h and .cpp files should be documented with doxygen comments (@file)". Affected files:'
            sh "grep '\\\\file\\|\\@file' -Lr . | grep '\\.cpp\\|\\.h'"
            sh "exit 1"
        }

        // check that no file contains NULL or assert
        try {
            // if .cpp or .h files contain NULL or assert, return 2
            sh "grep -lrE '([^_]NULL|[^_]assert)' . | grep -q '\\.cpp\\|\\.h' && exit 2 || exit 0"
        } catch (Exception e) {
            // change detected
            echo 'Usage of NULL and assert is prohibited. Affected files:'
            sh "grep -lrE '([^_]NULL|[^_]assert)' . | grep '\\.cpp\\|\\.h'"
            sh "exit 1"
        }

        // check that no file contains an #include <autopas/...>
        try {
            // if any file contains #include <autopas, return 2
            sh "grep -qlrE '#include <autopas' . && exit 2 || exit 0"
        } catch (Exception e) {
            // change detected
            echo 'Usage of #include <autopas...> is discouraged, please use "". Affected files:'
            sh "grep -lrE '#include <autopas' ."
            sh "exit 1"
        }

        // check that no file contains a typedef
        try {
            // if any file contains a typedef, return 2
            sh "grep -qlrE 'typedef ' . && exit 2 || exit 0"
        } catch (Exception e) {
            // change detected
            echo 'Usage of typedef is discouraged, please use using declarations instead. Affected files:'
            sh "grep -lrE 'typedef ' ."
            sh "exit 1"
        }
    }
}
