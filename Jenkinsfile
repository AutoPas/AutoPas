pipeline{
    agent {
        node {
            //cloud 'kubernetes'
            label 'openshift-autoscale'
        }
    }
    stages{
        stage('setup'){
            steps{
                echo 'Starting AutoPas Pipeline'
            }
        }
        stage("style check") {
            steps{
                parallel(
                    "build documentation": {
                        container('autopas-cmake-doxygen-make'){
                            dir("build-doxygen") {
                                sh 'cmake ..'
                                sh 'make doc_doxygen 2>DoxygenWarningLog.txt'
                            }
                        }
                        stash includes: 'build-doxygen/doc_doxygen/html/**', name: 'doxydocs'

                        // get doxygen warnings
                        recordIssues filters: [excludeFile('.*README.*')], tools: [doxygen(pattern: 'build-doxygen/DoxygenWarningLog.txt')], unstableTotalAll: 1
                    },
                    "clang format": {
                        dir("clang-format"){
                            container('autopas-clang6-cmake-ninja-make'){
                                sh "CC=clang CXX=clang++ cmake -G Ninja -DOPENMP=ON .."
                                sh "ninja clangformat"
                            }
                            script{
                                // return 2 if files have been modified by clang-format, 0 otherwise
                                try{
                                    // if files were modified, return 2
                                    sh "git status | grep -q modified && exit 2 || exit 0"
                                } catch (Exception e) {
                                    // change detected
                                    echo 'clang format errors detected. please format the code properly. Affected files:'
                                    sh "git status | grep modified"
                                    sh "exit 1"
                                }
                            }
                        }
                    },
                    "custom checks": {
                        dir("src"){
                            script{
                                // check if all header files have a #pragma once
                                try{
                                    // if header files do not contain #pragma once, make build unstable
                                    sh 'grep -L "#pragma once" -r . | grep -q "\\.h" && exit 2 || exit 0'
                                } catch (Exception e) {
                                    // change detected
                                    echo 'all header include guards should be implemented using "#pragma once". Affected files:'
                                    sh 'grep -L "#pragma once" -r . | grep "\\.h"'
                                    sh "exit 1"
                                }

                                // check if all files are documented with @file or \file doxygen comments
                                try{
                                    // if .cpp or .h files do not contain a file comment, return 2
                                    sh "grep '\\\\file\\|\\@file' -Lr . | grep -q '\\.cpp\\|\\.h' && exit 2 || exit 0"
                                } catch (Exception e) {
                                    // change detected
                                    echo 'all .h and .cpp files should be documented with doxygen comments (@file)". Affected files:'
                                    sh "grep '\\\\file\\|\\@file' -Lr . | grep '\\.cpp\\|\\.h'"
                                    sh "exit 1"
                                }
                            }
                        }
                    }
                )
            }
        }
        stage('build'){
            parallel{
                stage('gpu cloud') {
                    agent { label 'openshift-autoscale-gpu' }
                    steps{
                        container('cuda-10') {
                            dir("build"){
                                sh "cmake -DENABLE_CUDA=ON .."
                                sh "make -j 4 > buildlog-cuda.txt 2>&1 || (cat buildlog-cuda.txt && exit 1)"
                                sh "cat buildlog-cuda.txt"
                                sh "./tests/testAutopas/runTests"
                                sh "ctest -C checkExamples -j8 --verbose"
                            }
                        }
                    }
                    post{
                        always{
                            warnings canComputeNew: false, categoriesPattern: '', defaultEncoding: '', excludePattern: '', healthy: '', includePattern: '', messagesPattern: '', parserConfigurations: [[parserName: 'GNU Make + GNU C Compiler (gcc)', pattern: 'build*/buildlog-cuda.txt']], unHealthy: '', unstableTotalAll: '0', unstableTotalHigh: '0', unstableTotalLow: '0', unstableTotalNormal: '0'
                        }
                    }
                }
                stage("default") {
                    steps{
                        container('autopas-gcc7-cmake-make') {
                            dir("build"){
                                sh "cmake .."
                                sh "make -j 4 > buildlog.txt 2>&1 || (cat buildlog.txt && exit 1)"
                            }
                        }
                    }
                }
                stage("gcc openmp") {
                    steps{
                        container('autopas-gcc7-cmake-make') {
                            dir("build-openmp"){
                                sh "cmake -DOPENMP=ON .."
                                sh "make -j 4 > buildlog.txt 2>&1 || (cat buildlog.txt && exit 1)"
                            }
                        }
                    }
                }
                stage("gcc openmp address-sanitizer") {
                    steps{
                        container('autopas-gcc7-cmake-make') {
                            dir("build-openmp-address-sanitizer"){
                                sh "cmake -DOPENMP=ON -DCMAKE_BUILD_TYPE=Debug -DENABLE_ADDRESS_SANITIZER=ON .."
                                sh "make -j 4 > buildlog.txt 2>&1 || (cat buildlog.txt && exit 1)"
                            }
                        }
                    }
                }
                stage("address sanitizer") {
                    steps{
                        container('autopas-gcc7-cmake-make') {
                            dir("build-addresssanitizer"){
                                sh "cmake -DCMAKE_BUILD_TYPE=Debug -DENABLE_ADDRESS_SANITIZER=ON .."
                                sh "make -j 4 > buildlog.txt 2>&1 || (cat buildlog.txt && exit 1)"
                            }
                        }
                    }
                }
                stage("address sanitizer release") {
                    steps{
                        container('autopas-gcc7-cmake-make') {
                            dir("build-addresssanitizer-release"){
                                sh "cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_ADDRESS_SANITIZER=ON .."
                                sh "make -j 4 > buildlog.txt 2>&1 || (cat buildlog.txt && exit 1)"
                            }
                        }
                    }
                }
                stage("thread sanitizer") {
                    steps{
                        container('autopas-gcc7-cmake-make') {
                            dir("build-threadsanitizer"){
                                // this is for simple testing of our threading libraries.
                                sh "cmake -DCMAKE_BUILD_TYPE=Debug -DENABLE_THREAD_SANITIZER=ON .."
                                sh "make -j 4 > buildlog.txt 2>&1 || (cat buildlog.txt && exit 1)"
                            }
                        }
                    }
                }
                stage("clang openmp") {
                    steps{
                        container('autopas-clang6-cmake-ninja-make'){
                            dir("build-clang-ninja-openmp"){
                                sh "CC=clang CXX=clang++ cmake -G Ninja -DOPENMP=ON .."
                                sh "ninja -j 4 > buildlog_clang.txt 2>&1 || (cat buildlog_clang.txt && exit 1)"
                            }
                        }
                    }
                }
                stage("clang ninja address sanitizer") {
                    steps{
                        container('autopas-clang6-cmake-ninja-make'){
                            dir("build-clang-ninja-addresssanitizer-debug"){
                                sh "CXXFLAGS=-Wno-pass-failed CC=clang CXX=clang++ cmake -G Ninja -DCMAKE_BUILD_TYPE=Debug -DENABLE_ADDRESS_SANITIZER=ON .."
                                sh "ninja -j 4 > buildlog_clang.txt 2>&1 || (cat buildlog_clang.txt && exit 1)"
                            }
                        }
                    }
                }
                stage("clang ninja address sanitizer release") {
                    steps{
                        container('autopas-clang6-cmake-ninja-make'){
                            dir("build-clang-ninja-addresssanitizer-release"){
                                sh "CXXFLAGS=-Wno-pass-failed CC=clang CXX=clang++ cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -DENABLE_ADDRESS_SANITIZER=ON .."
                                sh "ninja -j 4 > buildlog_clang.txt 2>&1 || (cat buildlog_clang.txt && exit 1)"
                            }
                        }
                    }
                }
                stage("archer") {
                    steps{
                        container('autopas-archer'){
                            dir("build-archer"){
                                sh "CXXFLAGS=-Wno-pass-failed CC=clang-archer CXX=clang-archer++ cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -DUSE_VECTORIZATION=OFF .."
                                sh "ninja -j 4 > buildlog_clang.txt 2>&1 || (cat buildlog_clang.txt && exit 1)"
                            }
                        }
                    }
                }
                stage("intel") {
                    steps{
                        container('autopas-intel18'){
                            dir("build-intel"){
                                sh "bash -i -c 'which icc && CC=`which icc` CXX=`which icpc` cmake -DOPENMP=OFF ..'"
                                sh "bash -i -c 'make -j 8 > buildlog_intel.txt 2>&1 || (cat buildlog_intel.txt && exit 1)'"
                            }
                        }
                    }
                }
                stage("intel openmp") {
                    steps{
                        container('autopas-intel18'){
                            dir("build-intel-ninja-openmp"){
                                sh "bash -i -c 'which icc && CC=`which icc` CXX=`which icpc` cmake -G Ninja -DOPENMP=ON ..'"
                                sh "bash -i -c 'ninja -j 8 > buildlog_intel.txt 2>&1 || (cat buildlog_intel.txt && exit 1)'"
                            }
                        }
                    }
                }
            }
            post{
                always{
                    warnings canComputeNew: false, categoriesPattern: '', defaultEncoding: '', excludePattern: '', healthy: '', includePattern: '', messagesPattern: '', parserConfigurations: [[parserName: 'Clang (LLVM based)', pattern: 'build*/buildlog_clang.txt'], [parserName: 'GNU Make + GNU C Compiler (gcc)', pattern: 'build*/buildlog.txt'], [parserName: 'Intel C Compiler', pattern: 'build*/buildlog_intel.txt']], unHealthy: '', unstableTotalAll: '0', unstableTotalHigh: '0', unstableTotalLow: '0', unstableTotalNormal: '0'
                }
            }
        }
        stage("test") {
            steps{
                parallel(
                    "default": {
                        container('autopas-gcc7-cmake-make') {
                            dir("build"){
                                //sh "env CTEST_OUTPUT_ON_FAILURE=1 make test"
                                sh 'env GTEST_OUTPUT="xml:$(pwd)/test.xml" ./tests/testAutopas/runTests'
                            }
                        }
                    },
                    "gcc openmp": {
                        container('autopas-gcc7-cmake-make') {
                            dir("build-openmp"){
                                sh './tests/testAutopas/runTests'
                            }
                        }
                    },
                    "gcc openmp address-sanitizer": {
                        container('autopas-gcc7-cmake-make') {
                            dir("build-openmp-address-sanitizer"){
                                sh './tests/testAutopas/runTests'
                            }
                        }
                    },
                    "address sanitizer": {
                        container('autopas-gcc7-cmake-make') {
                            dir("build-addresssanitizer"){
                                sh './tests/testAutopas/runTests'
                            }
                        }
                    },
                    "address sanitizer release": {
                        container('autopas-gcc7-cmake-make') {
                            dir("build-addresssanitizer-release"){
                                sh './tests/testAutopas/runTests'
                            }
                        }
                    },
                    "thread sanitizer": {
                        container('autopas-gcc7-cmake-make') {
                            dir("build-threadsanitizer"){
                                sh './tests/testAutopas/runTests'
                            }
                        }
                    },
                    /*dir("build-memorysanitizer"){
                        sh './tests/testAutopas/runTests'
                    }*/
                    "clang openmp": {
                        container('autopas-clang6-cmake-ninja-make'){
                            dir("build-clang-ninja-openmp"){
                                sh './tests/testAutopas/runTests'
                            }
                        }
                    },
                    "clang ninja address sanitizer": {
                        container('autopas-clang6-cmake-ninja-make'){
                            dir("build-clang-ninja-addresssanitizer-debug"){
                                sh './tests/testAutopas/runTests'
                            }
                        }
                    },
                    "clang ninja address sanitizer release": {
                        container('autopas-clang6-cmake-ninja-make'){
                            dir("build-clang-ninja-addresssanitizer-release"){
                                sh './tests/testAutopas/runTests'
                            }
                        }
                    },
                    "archer": {
                        container('autopas-archer'){
                            dir("build-archer"){
                                // needed to properly check all test cases (also just single ones)
                                // I suspect that archer breaks when death tests are used
                                sh 'export TSAN_OPTIONS="ignore_noninstrumented_modules=1" && export ARCHER_OPTIONS="print_ompt_counters=1" && ctest --verbose'
                            }
                        }
                    },
                    "intel": {
                        container('autopas-intel18'){
                            dir("build-intel"){
                                sh "bash -i -c './tests/testAutopas/runTests'"
                            }
                        }
                    },
                    "intel openmp": {
                        container('autopas-intel18'){
                            dir("build-intel-ninja-openmp"){
                                sh "bash -i -c './tests/testAutopas/runTests'"
                            }
                        }
                    }
                )
            }
        }
        stage("checkExamples") {
            steps{
                parallel(
                    "default": {
                        container('autopas-gcc7-cmake-make') {
                            dir("build/examples") {
                                sh 'ctest -C checkExamples -j8 --verbose'
                            }
                        }
                    },
                    "gcc openmp": {
                        container('autopas-gcc7-cmake-make') {
                            dir("build-openmp/examples") {
                                sh 'ctest -C checkExamples -j8 --verbose'
                            }
                        }
                    },
                    /*"gcc openmp address-sanitizer": {
                        container('autopas-gcc7-cmake-make') {
                            dir("build-openmp-address-sanitizer/examples"){
                                sh 'ctest -C checkExamples -j8'
                            }
                        }
                    },*/
                    /*"address sanitizer": {
                        container('autopas-gcc7-cmake-make') {
                            dir("build-addresssanitizer/examples"){
                                sh 'ctest -C checkExamples -j8 --verbose'
                            }
                        }
                    },*/
                    /*"address sanitizer release": {
                        container('autopas-gcc7-cmake-make') {
                            dir("build-addresssanitizer-release/examples"){
                                sh 'ctest -C checkExamples -j8 --verbose'
                            }
                        }
                    },*/
                    /*"thread sanitizer": {
                        container('autopas-gcc7-cmake-make') {
                            dir("build-threadsanitizer/examples"){
                                sh 'ctest -C checkExamples -j8 --verbose'
                            }
                        }
                    },*/
                    "clang openmp": {
                        container('autopas-clang6-cmake-ninja-make'){
                            dir("build-clang-ninja-openmp/examples"){
                                sh 'ctest -C checkExamples -j8 --verbose'
                            }
                        }
                    },
                    /*"clang ninja address sanitizer": {
                        container('autopas-clang6-cmake-ninja-make'){
                            dir("build-clang-ninja-addresssanitizer-debug/examples"){
                                sh 'ctest -C checkExamples -j8 --verbose'
                            }
                        }
                    },*/
                    "clang ninja address sanitizer release": {
                        container('autopas-clang6-cmake-ninja-make'){
                            dir("build-clang-ninja-addresssanitizer-release/examples"){
                                sh 'ctest -C checkExamples -j8 --verbose'
                            }
                        }
                    },
                    "archer": {
                        container('autopas-archer'){
                            dir("build-archer/examples"){
                                sh 'export TSAN_OPTIONS="ignore_noninstrumented_modules=1" && export ARCHER_OPTIONS="print_ompt_counters=1" && ctest -C checkExamples -j8 --verbose'
                            }
                        }
                    },
                    "intel": {
                        container('autopas-intel18'){
                            dir("build-intel/examples"){
                                sh "bash -i -c 'ctest -C checkExamples -j8 --verbose'"
                            }
                        }
                    },
                    "intel openmp": {
                        container('autopas-intel18'){
                            dir("build-intel-ninja-openmp/examples"){
                                sh "bash -i -c 'ctest -C checkExamples -j8 --verbose'"
                            }
                        }
                    }
                )
            }
        }

        stage("generate reports"){
            steps{
                // get test results -- mainly to get number of tests
                junit 'build/test.xml'

                // generate coverage
                dir("coverage"){
                    container('autopas-build-code-coverage'){
                        sh "cmake -DCodeCoverage=ON -DCMAKE_BUILD_TYPE=Debug .."
                        sh "make AutoPas_cobertura -j 4"
                    }
                    cobertura autoUpdateHealth: false, autoUpdateStability: false, coberturaReportFile: 'coverage.xml', conditionalCoverageTargets: '70, 0, 0', failUnhealthy: false, failUnstable: false, lineCoverageTargets: '80, 0, 0', maxNumberOfBuilds: 0, methodCoverageTargets: '80, 0, 0', onlyStable: false, sourceEncoding: 'ASCII', zoomCoverageChart: false
                }
            }
        }

        stage("publish documentation"){
            when{ branch 'master' }
            agent{ label 'atsccs11_prio' }
            steps{
                unstash 'doxydocs'
                dir("build-doxygen"){
                    sh 'touch /import/www/wwwsccs/html/AutoPas/doxygen_doc/master || echo 0'
                    sh 'rm -rf /import/www/wwwsccs/html/AutoPas/doxygen_doc/master || echo 0'
                    sh 'cp -r doc_doxygen/html /import/www/wwwsccs/html/AutoPas/doxygen_doc/master'
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
