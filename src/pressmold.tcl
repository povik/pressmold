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
	{[-priority npriority_cuts] [-match_limit match_limit]}
proc prepare_cuts {args} {
	sta::parse_key_args "refine_mapping" args \
		keys {-priority -match_limit} \
		flags {}
	if {[info exists keys(-match_limit)]} {
		set match_limit $keys(-match_limit)
	} else {
		set match_limit 16
	}

	if {[info exists keys(-priority)]} {
		set priority $keys(-priority)
	} else {
		set priority 64
	}

	sta::prepare_cuts_cmd $priority $match_limit
}

sta::define_cmd_args "refine_mapping" \
	{[-seq pass_sequence] [-temp starting_temperature]}

proc refine_mapping {args} {
	sta::parse_key_args "refine_mapping" args \
		keys {-seq -temp} \
		flags {}

	if { ![info exists keys(-seq)] } {
		error "A -seq argument is required"
	}

	set seq $keys(-seq)
	set refs_blend 1.0
	for {set i 0} {$i < [string length $seq]} {incr i} {
		set crumb [string index $seq $i]

		set mark [expr $i + 1]
		while {($i < [string length $seq]) \
			&& [string is integer [string range $seq $mark [expr $i + 1]]]} {incr i}
		set rep [string range $seq $mark $i]
		if {![string length $rep]} {set rep 1}

		if {$crumb == "I"} {
			set round init
		} elseif {$crumb == "A"} {
			set round flow
		} elseif {$crumb == "E"} {
			set round exact
		} elseif {$crumb == "T"} {
			set round anneal
			if {![info exists keys(-temp)]} {
				error "Missing -temp argument for annealing"
			}
		} elseif {$crumb == "D"} {
			set round depth
		} elseif {$crumb == "d"} {
			set round depth2
		} elseif {$crumb == "S"} {
			set round save
		} elseif {$crumb == "s"} {
			set round stitch
		} else {
			error "Symbol $crumb passed in the -seq argument not recognized"
		}

		for {set j 0} {$j < $rep} {incr j} {
			if {$round == "anneal"} {
				if {$rep == 1} {
					set param $keys(-temp)
				} else {
					set param [expr (1.0 - ($j / ($rep - 1.0))) * $keys(-temp)]
				}
			} elseif {$round == "flow"} {
				set param $refs_blend
				set refs_blend [expr $refs_blend / 2.0]
			} else {
				set param [expr $i == 0]
			}
			sta::mapping_round_cmd $round $param
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
	refine_mapping -seq A5E3
	extract_mapping
}
