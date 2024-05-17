sta::define_cmd_args "register_all_cells" {[-verbose]}

proc register_cell {args} {
	sta::parse_key_args "register_cell" args \
		keys {} \
		flags {-verbose}

	sta::register_cell $args [info exists flags(-verbose)]
}

sta::define_cmd_args "register_all_cells" {[-verbose]}

proc register_all_cells {args} {
	sta::parse_key_args "register_all_cells" args \
		keys {} \
		flags {-verbose}

	foreach cell [get_lib_cells *] {
		sta::register_cell $cell [info exists flags(-verbose)]
	}
}

proc write_eqy_problem {path} {
	set top [[[sta::top_instance] cell] name]
	set eqy_conf [open $path.eqy w]

	puts $eqy_conf "\[gold]"
	sta::write_aig_verilog ${path}_gold.v $top
	puts $eqy_conf "read_verilog ${path}_gold.v\n"

	puts $eqy_conf "\[gate]"
	set lib_no 1
	foreach lib [get_libs *] {
		set fn ${path}_lib${lib_no}.lib
		sta::write_liberty $lib $fn
		puts $eqy_conf "read_liberty -ignore_miss_func $fn"
		set lib_no [expr $lib_no + 1]
	}
	write_verilog ${path}_mapping.v
	puts $eqy_conf "read_verilog ${path}_mapping.v"
	puts $eqy_conf "hierarchy -top $top\n"
	puts $eqy_conf "\[strategy basic]\nuse sat\ndepth 10\n"
	puts $eqy_conf "\[partition *]\namend *\n"
	close $eqy_conf
}

proc check_equivalence {path} {
	write_eqy_problem $path
	exec eqy -f ${path}.eqy >@stdout
}

sta::define_cmd_args "prepare_cuts" \
	{[-cuts cuts_limit] [-matches matches_limit] [-max_cut max_cut]}
proc prepare_cuts {args} {
	sta::parse_key_args "prepare_cuts" args \
		keys {-cuts -matches -max_cut} \
		flags {}

	if {[info exists keys(-matches)]} {
		set matches $keys(-matches)
	} else {
		set matches 16
	}

	if {[info exists keys(-cuts)]} {
		set cuts $keys(-cuts)
	} else {
		set cuts 64
	}

	if {[info exists keys(-max_cut)]} {
		set max_cut $keys(-max_cut)
	} else {
		# select the build-time limit
		set max_cut -1
	}

	sta::prepare_cuts_cmd $cuts $matches $max_cut
}

sta::define_cmd_args "develop_mapping" \
	{[-sequence pass_sequence] [-temperature starting_temperature]}

proc develop_mapping {args} {
	sta::parse_key_args "develop_mapping" args \
		keys {-sequence -temperature} \
		flags {}

	if { ![info exists keys(-sequence)] } {
		error "A -sequence argument is required"
	}

	set seq $keys(-sequence)
	set refs_blend 1.0
	for {set i 0} {$i < [string length $seq]} {incr i} {
		set crumb [string index $seq $i]

		set mark [expr $i + 1]
		while {($i < [string length $seq]) \
			&& [string is integer [string range $seq $mark [expr $i + 1]]]} {incr i}
		set rep [string range $seq $mark $i]
		if {![string length $rep]} {set rep 1}

		if {$crumb == "A"} {
			set round flow
		} elseif {$crumb == "E"} {
			set round exact
		} elseif {$crumb == "T"} {
			set round anneal
			if {![info exists keys(-temperature)]} {
				error "Missing -temperature argument for annealing"
			}
		} elseif {$crumb == "D"} {
			set round depth
		} elseif {$crumb == "d"} {
			set round depth2
		} elseif {$crumb == "X"} {
			set round fuzzy
		} elseif {$crumb == "S"} {
			set round save
		} elseif {$crumb == "s"} {
			set round stitch
		} else {
			error "Symbol $crumb passed in the -sequence argument not recognized"
		}

		set prev_round ""
		for {set j 0} {$j < $rep} {incr j} {
			set param2 false
			if {$round == "anneal"} {
				if {$rep == 1} {
					set param $keys(-temperature)
				} else {
					set param [expr (1.0 - ($j / ($rep - 1.0))) * $keys(-temperature)]
				}
			} elseif {$round == "flow"} {
				set param $refs_blend
				set refs_blend [expr $refs_blend / 2.0]
			} elseif {$round == "fuzzy"} {
				if {$rep == 1} {
					set param $keys(-temperature)
				} else {
					set param [expr (1.0 - ($j / ($rep - 0.0))) * $keys(-temperature)]
				}
				if {$prev_round != "fuzzy"} {
					set param2 true
				}
			} else {
				set param [expr $i == 0]
			}
			sta::mapping_round_cmd $round $param $param2
			set prev_round $round
		}
	}

	puts ""
}

proc read_aiger {path} {
	sta::read_aiger_cmd "$path" "top"
}

proc extract_mapping {} {
	sta::extract_mapping
}

proc report_aig {} {
	sta::report_aig
}

proc lose_choices {} {
	sta::lose_choices
}

proc report_mapping {} {
	sta::report_mapping
}

proc map {} {
	prepare_cuts
	develop_mapping -sequence A5E3
	extract_mapping
}

proc report_sibling_usage {} {
	sta::report_sibling_usage
}
