# Read in function that converts fasta to nexus
source("../code/seq_fxns.R")

# Read in functions that generate hard constraints and
# partial constraints from an input tree
source("../code/constraint_fxns.R")

# Generate nexus data file
# ========================
fasta2nexus("expanded_chalcidoidea1_aligned_trimmed.fasta","mb_runs/chalcidoidea1.nex")
fasta2nexus("expanded_chalcidoidea2_aligned_trimmed.fasta","mb_runs/chalcidoidea2.nex")
fasta2nexus("expanded_chalcidoidea3_aligned_trimmed.fasta","mb_runs/chalcidoidea3.nex")

# Generate chalcidoidea partial constraints
# =========================================
gen_mb_con_file("cruaud_CO1.nwk", "mb_runs/cruaud_constraints.nex")

# Generate higher constraints
# ===========================
out_file1 <- "mb_runs/higher_constraints1.nex"
out_file2 <- "mb_runs/higher_constraints2.nex"
out_file3 <- "mb_runs/higher_constraints3.nex"

# Read in taxonomy metadata
T1 <- read.delim("expanded_chalcidoidea1_taxonomy.tsv")
T2 <- read.delim("expanded_chalcidoidea2_taxonomy.tsv")
T3 <- read.delim("expanded_chalcidoidea3_taxonomy.tsv")

# Print header to output file
cat("#NEXUS\n\nbegin mrbayes;\n",file=out_file1)
cat("#NEXUS\n\nbegin mrbayes;\n",file=out_file2)
cat("#NEXUS\n\nbegin mrbayes;\n",file=out_file3)

# Assemble all tip labels and higher classification info
ingroup1 <- T1$TipLabel[T1$Family!="Mymarommatidae"]
ingroup2 <- T2$TipLabel[T2$Family!="Mymarommatidae"]
ingroup3 <- T3$TipLabel[T3$Family!="Mymarommatidae"]
add_hard_constraint("Chalcidoidea",ingroup1,out_file1)
add_hard_constraint("Chalcidoidea",ingroup2,out_file2)
add_hard_constraint("Chalcidoidea",ingroup3,out_file3)

# Print tail to output file
cat ("end;\n", file=out_file1, append=TRUE)
cat ("end;\n", file=out_file2, append=TRUE)
cat ("end;\n", file=out_file3, append=TRUE)

# Generate run files
# ==================
make_run_file <- function(ver) {

    out_file <- paste0("mb_runs/run",ver,".nex")
    
    cat("#NEXUS\n\nbegin mrbayes;\n",file=out_file)
    
    output <- function(a) { cat(a,"\n", sep="",file=out_file,append=TRUE) }
    
    output("\tset autoclose=yes nowarn=yes;")
    output('\tset dir="../";')
    exe_nex_file <- paste0("\texe expanded_chalcidoidea",ver,".nex;")
    output(exe_nex_file)
    output("\texe cruaud_constraints.nex;")
    exe_higher_file <- paste0("\texe higher_constraints",ver,".nex;")
    output(exe_higher_file)
    output('\tset dir="";')
    output("")
    output("\tcharset 1st = 1-. \3;")
    output("\tcharset 2nd = 2-. \3;")
    output("\tcharset 3rd = 3-. \3;")
    output("\tpartition cod = 3: 1st, 2nd, 3rd;")
    output("\tset partition=cod;")
    output("")
    output("\tprset brlenspr=clock:uniform clockvarpr=strict;")
    output("\tprset topologypr=constraints(1-505);")
    output("")
    output("\tlset rates=gamma nst=mixed;")
    output("\tprset ratepr=variable;")
    output("\tunlink shape=(all) statefreq=(all) revmat=(all);")
    output("")
    output("\tmcmc;")
    output("end;")
}

make_run_file(1)
make_run_file(2)
make_run_file(3)
