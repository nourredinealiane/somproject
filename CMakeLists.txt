cmake_minimum_required(VERSION 2.6)
project(som)

find_package(Qt5Widgets REQUIRED)
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
find_package(Threads REQUIRED)

# mettre ici les fichiers classés par types
set(SOURCES  windowSom.cpp sortingbox.cpp shapeitem.cpp math.cpp
	     vecteur.cpp  evaluation.cpp som.cpp qdata.cpp main.cpp)
set(HEADERS  windowSom.h sortingbox.h shapeitem.h math.h
             evaluation.h  vecteur.h  som.h qdata.h type.h)
#set(FORMS window.ui)

set(CMAKE_CXX_FLAGS ${CMAKE_C_FLAGS} "-std=c++0x ")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS}")

add_definitions(-Wall ${QT_DEFINITIONS} ${Qt5Widgets_DEFINITIONS})

include_directories(${Qt5Widgets_INCLUDE_DIRS} ${CMAKE_CURRENT_BINARY_DIR})
#include_directories(${Qt5Widgets_INCLUDE_DIRS} ${CMAKE_CURRENT_BINARY_DIR} -/usr/lib/libpq.so)

qt5_wrap_cpp(HEADERS_MOC ${HEADERS})

#qt5_wrap_ui(FORMS_MOC ${FORMS})

add_executable(som ${SOURCES} ${HEADERS_MOC} )#${FORMS_MOC})

# Use the Widgets module from Qt 5.

qt5_use_modules(som Sql Widgets)
