find -name "CMakeLists.txt" ! -path "./libs*" -exec cmake-format -i {} \;
