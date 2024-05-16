%include "tcl/Exception.i"
%include "tcl/StaTclTypes.i"
%include "tcl/StaTcl.i"
%include "tcl/NetworkEdit.i"
%include "sdf/Sdf.i"
%include "dcalc/DelayCalc.i"
%include "parasitics/Parasitics.i"
%include "power/Power.i"
%include "verilog/Verilog.i"

%{
	#include "commands.h"	
%}
extern bool register_cell_cmd(LibertyCell *cell, bool verbose);
extern void prepare_cuts_cmd(int npriority_cuts, int nmatches_max);
extern void read_aiger_cmd(const char *filename, const char *name);
extern void portlist_cmd();
extern void write_aig_verilog(const char *filename, const char *module_name);
extern void report_mapping();
extern void extract_mapping();
extern void mapping_round_cmd(const char *kind, float param, bool param2);
extern void report_aig();
extern void lose_choices();
extern void report_sibling_usage();
