pipeline{
    agent { label 'atsccs11' }
    stages{
        stage('setup'){
            steps{
                echo 'Starting AutoPas Pipeline'
                githubNotify context: 'build', description: 'build pending...',  status: 'PENDING', targetUrl: currentBuild.absoluteUrl
                githubNotify context: 'test', description: 'test pending...',  status: 'PENDING', targetUrl: currentBuild.absoluteUrl
            }
        }
        stage("build") {
            steps{
                githubNotify context: 'build', description: 'build in progress...',  status: 'PENDING', targetUrl: currentBuild.absoluteUrl
                dir("build"){
                    sh "cmake .."
                    sh "make"
                }
                dir("build-addresssanitizer"){
                    sh "cmake -DCodeCoverage=OFF -DCMAKE_BUILD_TYPE=Debug -DENABLE_ADDRESS_SANITIZER=ON .."
                    sh "make -j 4"
                }
                dir("build-threadsanitizer"){
                    // this is for simple testing of our threading libraries.
                    sh "cmake -DCodeCoverage=OFF -DCMAKE_BUILD_TYPE=Debug -DENABLE_THREAD_SANITIZER=ON .."
                    sh "make -j 4"
                }
            }
            post{
                success{
                    githubNotify context: 'build', description: currentBuild.durationString,  status: 'SUCCESS', targetUrl: currentBuild.absoluteUrl
                }
                failure{
                    githubNotify context: 'build', description: currentBuild.description, status: 'FAILURE', targetUrl: currentBuild.absoluteUrl
                }
                unstable{
                    githubNotify context: 'build', description: currentBuild.description, status: 'FAILURE', targetUrl: currentBuild.absoluteUrl
                }
                aborted{
                    githubNotify context: 'build', description: 'build aborted',  status: 'ERROR', targetUrl: currentBuild.absoluteUrl
                }
            }
        }
        stage("test") {
            steps{
                githubNotify context: 'test', description: 'test in progress...',  status: 'PENDING', targetUrl: currentBuild.absoluteUrl
                dir("build"){
                    //sh "env CTEST_OUTPUT_ON_FAILURE=1 make test"
                    sh 'env GTEST_OUTPUT="xml:$(pwd)/test.xml" ./tests/testAutopas/runTests'
                }
                dir("build-addresssanitizer"){
                    sh 'env GTEST_OUTPUT="xml:$(pwd)/test-address.xml" ./tests/testAutopas/runTests'
                }
                dir("build-threadsanitizer"){
                    sh 'env GTEST_OUTPUT="xml:$(pwd)/test-thread.xml" ./tests/testAutopas/runTests'
                }
            }
            post{
                success{
                    githubNotify context: 'test', description: currentBuild.durationString,  status: 'SUCCESS', targetUrl: currentBuild.absoluteUrl
                }
                failure{
                    githubNotify context: 'test', description: currentBuild.description,  status: 'FAILURE', targetUrl: currentBuild.absoluteUrl
                }
                unstable{
                    githubNotify context: 'test', description: currentBuild.description,  status: 'FAILURE', targetUrl: currentBuild.absoluteUrl
                }
                aborted{
                    githubNotify context: 'test', description: 'build aborted',  status: 'ERROR', targetUrl: currentBuild.absoluteUrl
                }
            }
        }
        stage("build documentation"){
            steps{
                dir("build") { sh 'make doc_doxygen 2>DoxygenWarningLog.txt' }
            }
        }
        stage("update documentation"){
            when{ branch 'master' }
            steps{
                dir("build"){
                    sh 'touch /import/www/wwwsccs/html/AutoPas/doxygen_doc/master || echo 0'
                    sh 'rm -rf /import/www/wwwsccs/html/AutoPas/doxygen_doc/master || echo 0'
                    sh 'cp -r doc_doxygen/html /import/www/wwwsccs/html/AutoPas/doxygen_doc/master'
                }
            }
        }
        stage("generate coverage report"){
            steps{
                junit 'build/test.xml'
                junit 'build-addresssanitizer/test-address.xml'
                junit 'build-threadsanitizer/test-thread.xml'
                warnings canComputeNew: false, canResolveRelativePaths: false, categoriesPattern: '', defaultEncoding: '', excludePattern: '.*README.*', healthy: '', includePattern: '', messagesPattern: '', parserConfigurations: [[parserName: 'Doxygen', pattern: 'build/DoxygenWarningLog.txt']], unHealthy: '', unstableTotalAll: '0'

                dir("coverage"){
                    sh "cmake -DCodeCoverage=ON -DCMAKE_BUILD_TYPE=Debug .."
                    sh "make AutoPas_cobertura"
                    cobertura autoUpdateHealth: false, autoUpdateStability: false, coberturaReportFile: 'coverage.xml', conditionalCoverageTargets: '70, 0, 0', failUnhealthy: false, failUnstable: false, lineCoverageTargets: '80, 0, 0', maxNumberOfBuilds: 0, methodCoverageTargets: '80, 0, 0', onlyStable: false, sourceEncoding: 'ASCII', zoomCoverageChart: false
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
        cleanup {
            when { environment ignoreCase: true, name: 'currentBuild.result', value: 'SUCCESS' }
            githubNotify context: 'continuous-integration/jenkins/branch', description: currentBuild.durationString,  status: 'SUCCESS', targetUrl: currentBuild.absoluteUrl
        }
    }
}
