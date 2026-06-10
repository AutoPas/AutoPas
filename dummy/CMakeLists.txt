project(pmt-dummy)

include(AddSensor)
add_sensor(SENSOR_NAME "dummy" HEADER Dummy.h SRC_FILES Dummy.cpp)
