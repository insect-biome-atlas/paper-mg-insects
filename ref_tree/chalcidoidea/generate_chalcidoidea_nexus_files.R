# Read in function that converts fasta to nexus
source("../code/seq_fxns.R")

# Read in functions that generate hard constraints and
# partial constraints from an input tree
source("../code/constraint_fxns.R")

# Generate nexus data file
# ========================
fasta2nexus("expanded_chalcidoidea1_aligned.fasta","mb_runs/chalcidoidea1.nex")
fasta2nexus("expanded_chalcidoidea2_aligned.fasta","mb_runs/chalcidoidea2.nex")

# Generate chalcidoidea partial constraints
# =========================================
gen_mb_con_file("cruaud_CO1.nwk", "mb_runs/cruaud_constraints.nex")
num_constraints1 <- 203 # Interior nodes in the Cruaud CO1-matched tree
num_constraints2 <- 203

# Generate higher constraints
# ===========================
out_file1 <- "mb_runs/higher_constraints1.nex"
out_file2 <- "mb_runs/higher_constraints2.nex"

# Read in taxonomy metadata
T1 <- read.delim("expanded_chalcidoidea1_taxonomy.tsv")
T2 <- read.delim("expanded_chalcidoidea2_taxonomy.tsv")

# Print header to output file
cat("#NEXUS\n\nbegin mrbayes;\n",file=out_file1)
cat("#NEXUS\n\nbegin mrbayes;\n",file=out_file2)

# Assemble all tip labels and higher classification info
ingroup1 <- T1$TipLabel[T1$Family!="Mymarommatidae"]
ingroup2 <- T2$TipLabel[T2$Family!="Mymarommatidae"]
add_hard_constraint("Chalcidoidea",ingroup1,out_file1)
add_hard_constraint("Chalcidoidea",ingroup2,out_file2)
num_constraints1 <- num_constraints1 + 1
num_constraints2 <- num_constraints2 + 1

# Add constraints for all clades
for (clade in unique(T2$Clade)) {
    ingroup <- T2$TipLabel[T2$Clade==clade]
    if (length(ingroup)>1) {
        add_hard_constraint(clade,ingroup,out_file2)
        num_constraints2 <- num_constraints2 + 1
    }
}

# Print tail to output file
cat ("end;\n", file=out_file1, append=TRUE)
cat ("end;\n", file=out_file2, append=TRUE)

# Generate run files
# ==================
make_run_file <- function(ver, num_constraints) {

    out_file <- paste0("mb_runs/run",ver,".nex")
    
    cat("#NEXUS\n\nbegin mrbayes;\n",file=out_file)
    
    output <- function(a) { cat(a,"\n", sep="",file=out_file,append=TRUE) }
    
    output("\tset autoclose=yes nowarn=yes;")
    output('\tset dir="../";')
    exe_nex_file <- paste0("\texe chalcidoidea",ver,".nex;")
    output(exe_nex_file)
    output("\texe cruaud_constraints.nex;")
    exe_higher_file <- paste0("\texe higher_constraints",ver,".nex;")
    output(exe_higher_file)
    output('\tset dir="";')
    output("")
    output("\tcharset 1st = 1-. \\3;")
    output("\tcharset 2nd = 2-. \\3;")
    output("\tcharset 3rd = 3-. \\3;")
    output("\tpartition cod = 3: 1st, 2nd, 3rd;")
    output("\tset partition=cod;")
    output("")
    output("\tprset brlenspr=clock:uniform clockvarpr=strict;")
    constraint_stmt <- paste0("\tprset topologypr=constraints(1-",num_constraints,");")
    output(constraint_stmt)
    output("")
    output("\tlset rates=gamma nst=mixed;")
    output("\tprset ratepr=variable;")
    output("\tunlink shape=(all) statefreq=(all) revmat=(all);")
    output("")
    output("\tmcmc;")
    output("end;")
}

make_run_file(1, num_constraints1)
make_run_file(2, num_constraints2)
