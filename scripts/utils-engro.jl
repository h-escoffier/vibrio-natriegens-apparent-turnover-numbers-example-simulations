using COBREXA, DataFrames
import CSV
import AbstractFBCModels as A


function build_reaction_isozymes(
    model,
    wildkcat_path::String,
    proteome_path::String;
    scale::Real = 3.6 * 1e-3,
)
    
    clamp_kcat(x) = x === nothing ? nothing : clamp(x, 0.1, 1e7) # ??? 

    # 1. Load WkC file 
    df = CSV.read(wildkcat_path, DataFrame; delim = '\t')
    # rm missing values
    filter!(row -> !ismissing(row.rxn) && !ismissing(row.genes) &&
                   !ismissing(row.direction) && !ismissing(row.kcat), df)

    kcat_data = Dict{String,Dict{String,Tuple{Union{Float64,Nothing},Union{Float64,Nothing}}}}()
    for g in groupby(df, [:rxn, :genes])
        rxn  = g.rxn[1]
        gene = join(sort(split(g.genes[1], ';')), ';')
        fwd_vec = g.kcat[g.direction .== "forward"]
        rev_vec = g.kcat[g.direction .== "reverse"]

        get!(kcat_data, rxn, Dict())[gene] = (
            isempty(fwd_vec) ? nothing : clamp_kcat(fwd_vec[1]),
            isempty(rev_vec) ? nothing : clamp_kcat(rev_vec[1]),
        )
    end

    # load proteome data to get the stochiometry 
    proteome_data = Dict(String(row.gene) => (Float64(row.mw),
            ismissing(row.subunit_structure) ? "" : String(row.subunit_structure))
            for row in CSV.File(proteome_path; delim = '\t')
                if !ismissing(row.gene) && !ismissing(row.mw)
    )

    # create reaction isozyme object (based on supplementary materials of COBREXA2)
    mer_map = Dict("Homomonomer"=>1.0,"Monomer"=>1.0,"Homodimer"=>2.0,"Homotrimer"=>3.0,
                   "Homotetramer"=>4.0,"Homopentamer"=>5.0,"Homohexamer"=>6.0,
                   "Homoheptamer"=>7.0,"Homooctamer"=>8.0,"Homodecamer"=>10.0,"Homododecamer"=>12.0)

    reaction_isozymes = Dict{String,Dict{String,Isozyme}}()

    for rid in A.reactions(model)
        haskey(kcat_data, rid)       || continue
        has_reaction_grr(model, rid) || continue
        haskey(reaction_isozymes, rid) && continue # Keep only one isozyme

        grr = first(A.reaction_gene_association_dnf(model, rid))  # Keep first
        key       = join(sort(grr), ";")
        fwd, rev  = get(kcat_data[rid], key, get(kcat_data[rid], "kcat", (nothing, nothing)))

        ss_counts = if length(grr) == 1
            gid = first(grr)
            ss  = (proteome_data !== nothing && haskey(proteome_data, gid)) ?
                    get(mer_map, proteome_data[gid][2], 1.0) : 1.0
            [ss]
        else
            fill(1.0, length(grr))
        end

        reaction_isozymes[rid] = Dict("isozyme" => Isozyme(
            gene_product_stoichiometry = Dict(grr .=> ss_counts),
            kcat_forward = fwd === nothing ? nothing : fwd * scale, # Convert from s-1 to h-1
            kcat_reverse = rev === nothing ? nothing : rev * scale,
        ))
    end
    return reaction_isozymes
end
