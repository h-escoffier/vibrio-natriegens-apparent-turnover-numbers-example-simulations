using CSV, DataFrames, DataFramesMeta, Serialization, Statistics
using COBREXA, SBMLFBCModels, AbstractFBCModels

include("utils-engro.jl")

isdir("tmp") || mkdir("tmp") # create dir if necessary

# save molar mass data
mw = Dict(
    String(row.gene) => (Float64(row.mw))
    for row in CSV.File("data/engro2-molar-masses.tsv"; delim = '\t')
        if !ismissing(row.gene) && !ismissing(row.mw)
            )

serialize("tmp/molar-masses.js", mw)

# Unlike vnat analysis, the model is not condition specific
model = convert(
    AbstractFBCModels.CanonicalModel.Model, 
    AbstractFBCModels.load("data/model-engro2.xml")) 

strip_prefix(id, prefix) = startswith(id, prefix) ? id[length(prefix)+1:end] : id
model = identifier_map(
    model,
    reaction_map = id -> strip_prefix(id, "R_"),
    metabolite_map = id -> strip_prefix(id, "M_"),
    gene_map = id -> strip_prefix(id, "G_"),
    compartment_map = id -> strip_prefix(id, "C_"),
    coupling_map = id -> strip_prefix(id, "X_"),
);

serialize("tmp/model.js", model)

# kinetics paramaters (retrived using WILDkCAT) [units: M/h]
reaction_isozymes = build_reaction_isozymes(model, "data/kcat.tsv", "data/engro2-molar-masses.tsv")
serialize("tmp/reaction-isozymes.js", reaction_isozymes)

# currently the test is done only on mcf7 condition 

condition = "MCF7"

# NOT USED - serialize a specific flux sample [units: mmol/gDW/h]
# fluxes = Dict(
#     Symbol(r.Reaction) => r.Flux for r in CSV.File("data/"*condition*"-fluxes.csv")
# )
# serialize("tmp/"*condition*"-fluxes.js", fluxes)

# serialize the proteomics mass fraction data
table = CSV.File(
    "data/"*condition*"-mass-fraction.tsv";
    delim='\t'
    )

value_cols = filter(name -> name ∉ (:ENSEMBLID, :UniProtIDs), propertynames(table)) # Select columns 
mass = Dict(
    row.ENSEMBLID => mean(skipmissing([getproperty(row, col) for col in value_cols])) # or median() ?
    for row in table)

mass = Dict(
    String15(k) => Float64(v)
    for (k, v) in mass
    if v != 0
)

serialize("tmp/"*condition*"-mass-fractions.js", mass)

# serialize the proteomics absolute mole fraction data [units: ng/gDW]
table = CSV.File(
    "data/"*condition*"-absolute-mole-fraction.tsv";
    delim='\t'
    )

value_cols = filter(name -> name ∉ (:ENSEMBLID, :UniProtIDs), propertynames(table))
mole = Dict(
    row.ENSEMBLID => mean(skipmissing([getproperty(row, col) for col in value_cols])) # or median() ?
    for row in table)

mole = Dict(
    Symbol(k) => Float64(v)
    for (k, v) in mole
    if v != 0
)

serialize("tmp/"*condition*"-mole-fractions.js", mole)
