pipeline{
    agent { label 'atsccs11' }
    stages{
        stage('setup'){
            steps{
                echo 'Starting AutoPas Pipeline'
                githubNotify context: 'build', description: 'build pending...',  status: 'PENDING'
                githubNotify context: 'test', description: 'test pending...',  status: 'PENDING'
            }
        }
        stage("build") {
            steps{
                githubNotify context: 'build', description: 'build in progress...',  status: 'PENDING'
                dir("build"){
                    sh "cmake .."
                    sh "make"
                }
            }
            post{
                success{
                    githubNotify context: 'build', description: 'build successful! (${currentBuild.durationString})',  status: 'SUCCESS'
                }
                failure{
                    githubNotify context: 'build', description: 'build failed: ${currentBuild.description}',  status: 'FAILURE'
                }
                unstable{
                    githubNotify context: 'build', description: 'build unstable: ${currentBuild.description}',  status: 'FAILURE'
                }
                aborted{
                    githubNotify context: 'build', description: 'build aborted',  status: 'ERROR'
                }
            }
        }
        stage("test") {
            steps{
                githubNotify context: 'test', description: 'test in progress...',  status: 'PENDING'
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
            post{
                success{
                    githubNotify context: 'test', description: 'test successful! (${currentBuild.durationString})',  status: 'SUCCESS'
                }
                failure{
                    githubNotify context: 'test', description: 'test failed: ${currentBuild.description}',  status: 'FAILURE'
                }
                unstable{
                    githubNotify context: 'test', description: 'test unstable: ${currentBuild.description}',  status: 'FAILURE'
                }
                aborted{
                    githubNotify context: 'test', description: 'build aborted',  status: 'ERROR'
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
