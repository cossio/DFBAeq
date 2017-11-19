# DFBAeq
Scripts used to simulate steady states of genome-scale metabolic networks in continuous cell culture.
See http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005835.

# Data files
 - `IMDM.txt`. Concentrations of Iscove modified Dulbecco's medium.
 - `iCHOv1_K1_final.xml`. iCHO-K1 metabolic network. Available as supporting information in this paper: https://www.ncbi.nlm.nih.gov/pubmed/27883890.
 - `Shlomi.txt`. Contains kcat and molecular weights of enzymes in the network. Original source: http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002018.
 
# Scripts
 -  `c.py`. Generates the medium definition file, `c.txt`.
 - `cell.py`. Generates the metabolic network in the format that we use here, that is, the files `cho.mets`, `cho.rxns` and `cho.sto`.
 - `eq.jl`. Solves the optimization problem.

# Usage

The inputs to the script `eq.jl` are:
 - `cho.mets`, where each row is a metabolite, and the columns are: `id` (metabolite id), `y` (biomass coefficient), `e` (maintenance demand), `L` (secretion bound), `V` (uptake bound)`.
 - `cho.rxns`, where each row is a reaction, and the columns are: `id` (reaction id), `lb,ub` (flux lower and upper bounds), `ap,an` (flux enzymatic costs, in forward and backward directions).
 - `cho.sto`, stoichiometric matrix, where each row is a metabolite.

# Requirements
Python packages:
 - COBRApy (https://github.com/opencobra/cobrapy)
 - scipy, argparse, pandas (available through `pip`)

Julia packages:
 - DataFrames.jl, Gurobi.jl, ArgParse.jl (install with `Pkg.add("PkgName")`)


If you have any questions or have trouble using the scripts, please file an issue.
