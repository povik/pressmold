# Press Mold Mapper: OpenSTA-based standard cell mapper

Press Mold Mapper is an experimental standard cell mapper built on top of OpenSTA and derived from [toymap](https://github.com/povik/toymap).

## Sample session

### Library preparation

	read_liberty sky130_fd_sc_hd__tt_025C_1v80.lib
	register_all_cells

Alternatively you can cherry pick, e.g.

	foreach cell [get_lib_cells *_0] {
		register_cell -verbose $cell
	}
	foreach cell [get_lib_cells *_1] {
		register_cell -verbose $cell
	}

to only register those cells ending in `_0` and `_1`.

### Reading the design

	read_aiger subjects/priority_2.aig
	report_aig

outputs

```
section 'n' (110): ignoring
Read network 'top' with 1159 nodes
Mapping problem summary:
  128 inputs 8 outputs 1159 nodes 147 choice pairs 265 xor/mux detections
```

### Mapping

To ask for 5 rounds of area flow, 3 rounds of exact area passes:

	prepare_cuts
	refine_mapping -seq A5E3
	report_mapping

outputs

```
Cut matching statistics:
  1159 nodes 0.01 MiB cut cache 1.31 MiB match cache
  saturated 839 cuts (72.4 %), 258 matches (22.3 %)
  matches 10.3 mean 11.1 geom

  flow  A=  3219.3  N=  665  E= 1714  (blend=1.000)
  flow  A=  2894.0  N=  614  E= 1577  (blend=0.500)
  flow  A=  2710.1  N=  574  E= 1489  (blend=0.250)
  flow  A=  2580.0  N=  550  E= 1428  (blend=0.125)
  flow  A=  2554.9  N=  541  E= 1412  (blend=0.062)
 exact  A=  2419.8  N=  463  E= 1226
 exact  A=  2229.6  N=  396  E= 1021
 exact  A=  2067.0  N=  331  E=  891


   sky130_fd_sc_hd__a211oi_1     2  1.501e+01
   sky130_fd_sc_hd__a21boi_0     1  7.507e+00
     sky130_fd_sc_hd__a21o_1     1  7.507e+00
    sky130_fd_sc_hd__a21oi_1    19  9.509e+01
   sky130_fd_sc_hd__a221oi_1     1  8.758e+00
    sky130_fd_sc_hd__a31oi_1     3  1.877e+01
     sky130_fd_sc_hd__and2_0     2  1.251e+01
     sky130_fd_sc_hd__and3_1     3  1.877e+01
  sky130_fd_sc_hd__lpflow_isobufsrc_1     5  3.128e+01
     sky130_fd_sc_hd__mux2_1    67  7.545e+02
    sky130_fd_sc_hd__mux2i_1     3  3.003e+01
    sky130_fd_sc_hd__nand2_1    43  1.614e+02
   sky130_fd_sc_hd__nand2b_1     5  3.128e+01
    sky130_fd_sc_hd__nand3_1     4  2.002e+01
    sky130_fd_sc_hd__nand4_1     1  6.256e+00
   sky130_fd_sc_hd__nand4b_1     1  8.758e+00
     sky130_fd_sc_hd__nor2_1    83  3.115e+02
     sky130_fd_sc_hd__nor3_1    26  1.301e+02
     sky130_fd_sc_hd__nor4_1    14  8.758e+01
  sky130_fd_sc_hd__o2111ai_1     1  8.758e+00
   sky130_fd_sc_hd__o211ai_1     1  7.507e+00
    sky130_fd_sc_hd__o21ai_0    29  1.451e+02
   sky130_fd_sc_hd__o21bai_1     1  7.507e+00
   sky130_fd_sc_hd__o221ai_1     1  8.758e+00
    sky130_fd_sc_hd__o22ai_1     1  6.256e+00
     sky130_fd_sc_hd__o31a_1     1  8.758e+00
    sky130_fd_sc_hd__o31ai_1     3  2.252e+01
    sky130_fd_sc_hd__o32ai_1     1  8.758e+00
      sky130_fd_sc_hd__or2_0     2  1.251e+01
      sky130_fd_sc_hd__or3_1     2  1.251e+01
     sky130_fd_sc_hd__or3b_1     1  8.758e+00
      sky130_fd_sc_hd__or4_1     3  2.252e+01

Sum: 331 cells 2.037e+03 area
```

### Timing

	extract_mapping
	report_checks -unconstrained

outputs

```
Startpoint: A[119] (input port)
Endpoint: P[0] (output port)
Path Group: unconstrained
Path Type: max

  Delay    Time   Description
---------------------------------------------------------
   0.00    0.00 ^ input external delay
   0.00    0.00 ^ A[119] (in)
   0.02    0.02 v _00000239_/Y (sky130_fd_sc_hd__inv_1)
   0.18    0.20 ^ _00000260_/Y (sky130_fd_sc_hd__o21ai_0)
   0.09    0.29 v _00000261_/Y (sky130_fd_sc_hd__a21oi_1)
```
...
```
   0.17   20.40 ^ _00000336_/X (sky130_fd_sc_hd__mux2_1)
   0.07   20.48 v _00000338_/Y (sky130_fd_sc_hd__mux2i_1)
   0.19   20.67 ^ _00000339_/Y (sky130_fd_sc_hd__o21ai_0)
   0.10   20.77 v _00000342_/Y (sky130_fd_sc_hd__a21oi_1)
   0.19   20.96 ^ _00000343_/Y (sky130_fd_sc_hd__o21ai_0)
   0.10   21.05 v _00000346_/Y (sky130_fd_sc_hd__a21oi_1)
   0.18   21.24 ^ _00000347_/Y (sky130_fd_sc_hd__o21ai_0)
   0.07   21.30 v _00000348_/Y (sky130_fd_sc_hd__nand2_1)
   0.08   21.38 ^ _00000349_/Y (sky130_fd_sc_hd__o21ai_0)
   0.11   21.49 v _00000350_/Y (sky130_fd_sc_hd__mux2i_1)
   0.11   21.60 ^ _00000385_/X (sky130_fd_sc_hd__lpflow_isobufsrc_1)
   0.05   21.65 v _00000387_/Y (sky130_fd_sc_hd__nor3_1)
   0.07   21.72 ^ _00000388_/Y (sky130_fd_sc_hd__a21oi_1)
   0.00   21.72 ^ P[0] (out)
          21.72   data arrival time
---------------------------------------------------------
(Path is unconstrained)


```

### Verification

To check mapping equivalence with EQY (written out files will be prefixed with `/tmp/mapping_lec`):

	extract_mapping
	check_equivalence /tmp/mapping_lec

outputs

```
EQY  0:54:53 [mapping_lec] read_gold: starting process "yosys -ql mapping_lec/gold.log mapping_lec/gold.ys"
EQY  0:54:53 [mapping_lec] read_gold: finished (returncode=0)
EQY  0:54:53 [mapping_lec] read_gate: starting process "yosys -ql mapping_lec/gate.log mapping_lec/gate.ys"
EQY  0:54:53 [mapping_lec] read_gate: finished (returncode=0)
EQY  0:54:53 [mapping_lec] combine: starting process "yosys -ql mapping_lec/combine.log mapping_lec/combine.ys"
EQY  0:54:53 [mapping_lec] combine: finished (returncode=0)
EQY  0:54:53 [mapping_lec] partition: starting process "cd mapping_lec; yosys -ql partition.log partition.ys"
EQY  0:54:54 [mapping_lec] partition: finished (returncode=0)
```
...
```
EQY  0:54:55 [mapping_lec] Successfully proved designs equivalent
EQY  0:54:55 [mapping_lec] summary: Elapsed clock time [H:MM:SS (secs)]: 0:00:02 (2)
EQY  0:54:55 [mapping_lec] summary: Elapsed process time [H:MM:SS (secs)]: 0:00:02 (2)
EQY  0:54:55 [mapping_lec] DONE (PASS, rc=0)
```

### Writing the netlist

	extract_mapping
	write_verilog netlist.v


### Removing structural choices

To remove structural choices (will reset current mapping):

	lose_choices

outputs

	Cleared 147 choice pairs
	Removed 578 unused nodes

## Copyright notice

Except for the `third_party/OpenSTA` submodule, and the `src/cmake/swig_lib.cmake` file, the code is:

Copyright 2024 Martin Povišer

No explicit license assigned at this point
