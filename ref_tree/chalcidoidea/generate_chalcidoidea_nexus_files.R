# Read in function that converts fasta to nexus
source("../code/seq_fxns.R")

# Read in functions that generate hard constraints and
# partial constraints from an input tree
source("../code/constraint_fxns.R")

# Generate nexus data file
# ========================
fasta2nexus("expanded_chalcidoidea_aligned.fasta","mb_runs/chalcidoidea.nex")

# Generate chalcidoidea partial constraints
# =========================================
gen_mb_con_file("cruaud_CO1.nwk", "mb_runs/cruaud_constraints.nex")
num_constraints <- 203 # Interior nodes in the Cruaud CO1-matched tree

# Generate higher constraints
# ===========================
out_file <- "mb_runs/higher_constraints.nex"

# Read in taxonomy metadata
T <- read.delim("expanded_chalcidoidea_taxonomy.tsv")

# Print header to output file
cat("#NEXUS\n\nbegin mrbayes;\n",file=out_file)

# Assemble all tip labels and higher classification info
ingroup <- T$TipLabel[T$Family!="Mymarommatidae"]
add_hard_constraint("Chalcidoidea",ingroup,out_file)
num_constraints <- num_constraints + 1

# Add constraints for all clades
for (clade in unique(T$Clade)) {
    ingroup <- T$TipLabel[T$Clade==clade]
    if (length(ingroup)>1) {
        add_hard_constraint(clade,ingroup,out_file)
        num_constraints <- num_constraints + 1
    }
}

# Print tail to output file
cat ("end;\n", file=out_file, append=TRUE)

# Generate run file
# =================
make_run_file <- function(num_constraints) {

    out_file <- paste0("mb_runs/run.nex")
    
    cat("#NEXUS\n\nbegin mrbayes;\n",file=out_file)
    
    output <- function(a) { cat(a,"\n", sep="",file=out_file,append=TRUE) }
    
    output("\tset autoclose=yes nowarn=yes;")
    output('\tset dir="../";')
    output("\texe chalcidoidea.nex;")
    output("\texe cruaud_constraints.nex;")
    output("\texe higher_constraints.nex;")
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

make_run_file(num_constraints)
