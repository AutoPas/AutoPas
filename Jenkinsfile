pipeline{
    agent { label 'atsccs11' }
    stages{
        stage('setup'){
            steps{
                echo 'Starting AutoPas Pipeline'
            }
        }
        stage("build") {
            steps{
                dir("build"){
                    sh "cmake .."
                    sh "make"
                }
            }
        }
        stage("test") {
            steps{
                dir("build"){
                    //sh "env CTEST_OUTPUT_ON_FAILURE=1 make test"
                    sh 'env GTEST_OUTPUT="xml:$(pwd)/test.xml" ctest'
                }
                dir("build-addresssanitizer"){
                    sh "cmake -DCodeCoverage=OFF -DCMAKE_BUILD_TYPE=Debug -DENABLE_ADDRESS_SANITIZER=ON .."
                    sh "make -j 4"
                    sh "make test"
                }
                dir("build-threadsanitizer"){
                    sh "cmake -DCodeCoverage=OFF -DCMAKE_BUILD_TYPE=Debug -DENABLE_THREAD_SANITIZER=ON .."
                    sh "make -j 4"
                    sh "make test"
                }
            }
        }
        stage("documentation and coverage"){
            steps{
                dir("build"){
                    sh 'make doc_doxygen 2>DoxygenWarningLog.txt'
                    sh 'touch /import/www/wwwsccs/html/AutoPas/doxygen_doc/master || echo 0'
                    sh 'rm -rf /import/www/wwwsccs/html/AutoPas/doxygen_doc/master || echo 0'
                    sh 'cp -r doc_doxygen/html /import/www/wwwsccs/html/AutoPas/doxygen_doc/master'
                }
                junit 'build/test.xml'
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
    }
}
