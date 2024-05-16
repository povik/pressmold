set (CMAKE_CXX_STANDARD 20)
cmake_policy(SET CMP0078 NEW)
cmake_policy(SET CMP0086 NEW)

find_package(SWIG 4.0 REQUIRED)
include(UseSWIG)

swig_lib(NAME pressmold_swig
	NAMESPACE sta # TODO
	I_FILE pressmold.i
	SWIG_INCLUDES ${OPENSTA_HOME}
	SCRIPTS pressmold.tcl
		${OPENSTA_HOME}/tcl/Init.tcl
		${OPENSTA_HOME}/tcl/Util.tcl
		${OPENSTA_HOME}/tcl/CmdArgs.tcl
		${OPENSTA_HOME}/tcl/CmdUtil.tcl
		${OPENSTA_HOME}/tcl/Graph.tcl
		${OPENSTA_HOME}/tcl/Liberty.tcl
		${OPENSTA_HOME}/tcl/Link.tcl
		${OPENSTA_HOME}/tcl/Network.tcl
		${OPENSTA_HOME}/tcl/NetworkEdit.tcl
		${OPENSTA_HOME}/tcl/Property.tcl
		${OPENSTA_HOME}/tcl/Sdc.tcl
		${OPENSTA_HOME}/tcl/Search.tcl
		${OPENSTA_HOME}/tcl/Sta.tcl
		${OPENSTA_HOME}/tcl/Splash.tcl
		${OPENSTA_HOME}/tcl/Variables.tcl
		${OPENSTA_HOME}/tcl/WritePathSpice.tcl
		${OPENSTA_HOME}/verilog/Verilog.tcl
		${OPENSTA_HOME}/dcalc/DelayCalc.tcl
		${OPENSTA_HOME}/parasitics/Parasitics.tcl
		${OPENSTA_HOME}/power/Power.tcl
		${OPENSTA_HOME}/sdf/Sdf.tcl)

target_include_directories(pressmold_swig
	PUBLIC ${OPENSTA_HOME} ${OPENSTA_HOME}/include/sta/ .)

find_library(TCL_READLINE_LIBRARY tclreadline REQUIRED)

add_executable(pressmold
	pressmold.cc
	npn.cc)

target_include_directories(pressmold
	PUBLIC .)

target_link_libraries(pressmold
	PRIVATE pressmold_swig OpenSTA tclreadline)
