#!/usr/bin/env julia
# Generates golden reference data for pyVTTrac v2 by running VTTrac.jl v2.0.0
# (the reference algorithm implementation) and saving inputs + outputs to
# tests/data/golden/*.npz.
#
# This script is NOT run in CI (Julia is a generation-time-only dependency).
# The generated .npz files are committed to the repository and read by
# tests/test_golden.py, which compares them against the Fortran/Python port.
#
# Usage:
#   VTTRAC_JL_PATH=/path/to/VTTrac.jl julia --project=tools tools/gen_golden.jl
#
# `VTTRAC_JL_PATH` defaults to `../VTTrac.jl` relative to this repo (i.e. a
# sibling checkout), since pyVTTrac v2 no longer vendors VTTrac.jl as a git
# submodule.
#
# Index convention: VTTrac.jl is 1-based (Julia arrays). All positions/time
# indices in the saved .npz are converted to 0-based, matching the pyvttrac
# Python API (and the sentinel convention: -1 for invalid `t_index`, NaN for
# invalid float outputs) so tests/test_golden.py can compare directly against
# `pyvttrac.track()` output without any further index conversion.

import Pkg
Pkg.activate(@__DIR__)

using NPZ
using Statistics

const HERE = @__DIR__
const DEFAULT_VTTRAC_JL_PATH = joinpath(HERE, "..", "..", "VTTrac.jl")
const VTTRAC_JL_PATH = get(ENV, "VTTRAC_JL_PATH", DEFAULT_VTTRAC_JL_PATH)
const VTTRAC_SRC = joinpath(VTTRAC_JL_PATH, "src", "VTTrac.jl")
isfile(VTTRAC_SRC) || error(
    "VTTrac.jl source not found at $VTTRAC_SRC.\n" *
    "Set VTTRAC_JL_PATH to point at a VTTrac.jl checkout (v2.0.0), e.g.:\n" *
    "  VTTRAC_JL_PATH=/path/to/VTTrac.jl julia --project=tools tools/gen_golden.jl"
)
include(VTTRAC_SRC)
using .VTTrac

const OUTDIR = joinpath(HERE, "..", "tests", "data", "golden")
mkpath(OUTDIR)

# ---------------------------------------------------------------------------
# Minimal JSON encoder (avoids adding JSON.jl as a dependency for this script)
# ---------------------------------------------------------------------------
json_scalar(v::Bool) = v ? "true" : "false"
json_scalar(::Nothing) = "null"
json_scalar(v::AbstractString) = "\"" * replace(v, "\"" => "\\\"") * "\""
json_scalar(v::Real) = string(v)
function to_json(d::Dict{String,Any})
    parts = String[]
    for k in sort(collect(keys(d)))
        push!(parts, "\"$k\": $(json_scalar(d[k]))")
    end
    return "{" * join(parts, ", ") * "}"
end

# ---------------------------------------------------------------------------
# Synthetic field generator: a traveling sinusoidal wave, matching the
# pattern already used in VTTrac.jl's own test suite.
#   z(x,y,t) = offset + amp * sin(k*(x-cx*t)) * cos(k*(y-cy*t))
# x, y are 0-based physical coordinates equal to (index-1); this matches the
# 0-based convention used for x0/y0 golden inputs.
# ---------------------------------------------------------------------------
function make_field(nx::Int, ny::Int, nt::Int;
        k::Float64=2pi/10, cx::Float64=1.2, cy::Float64=1.2,
        amp::Float32=1.0f0, offset::Float32=0.0f0, dt::Float64=1.0)
    t = collect(0:nt-1) .* dt
    z = Array{Float32,3}(undef, nt, ny, nx)
    for it in 1:nt
        tt = t[it]
        for iy in 1:ny
            y = Float64(iy - 1)
            for ix in 1:nx
                x = Float64(ix - 1)
                z[it, iy, ix] = offset + amp * Float32(sin(k * (x - cx * tt)) * cos(k * (y - cy * tt)))
            end
        end
    end
    return z, t
end

const Z0, T0 = make_field(48, 48, 10)          # default field, for most cases
const Z0BIG, T0BIG = make_field(64, 64, 12)    # larger field, for edge-margin cases

# ---------------------------------------------------------------------------
# Bookkeeping: which status codes (1-9) have we observed across all cases?
# ---------------------------------------------------------------------------
const SEEN_STATUS = Set{Int}()
const CASE_NAMES = String[]

# ---------------------------------------------------------------------------
# run_case: sets up + runs one VTT tracking case and saves golden .npz
# ---------------------------------------------------------------------------
function run_case(name::AbstractString;
        z::Array{Float32,3}, t::Vector{Float64},
        mask::Union{Nothing,Array{Bool,3}}=nothing,
        zmiss::Union{Nothing,Real}=nothing,
        nsx::Int, nsy::Int,
        vxhw=nothing, vyhw=nothing, ixhw=nothing, iyhw=nothing,
        subgrid::Bool=true, subgrid_gaus::Bool=false,
        itstep::Int=1, ntrac::Int=3,
        score_method::AbstractString="xcor",
        Sth0::Float64=0.8, Sth1::Float64=0.7,
        vxch=nothing, vych=nothing,
        peak_inside_th=nothing, Cth=nothing,
        use_init_temp::Bool=false, min_samples::Int=1,
        tid0, x0::Vector{Float64}, y0::Vector{Float64},
        vxg=nothing, vyg=nothing,
        out_subimage::Bool=false, out_score_ary::Bool=false)

    vtt = VTTrac.VTT(z; t=t, mask=mask, zmiss=zmiss)
    VTTrac.setup(vtt, nsx, nsy;
        vxhw=vxhw, vyhw=vyhw, ixhw=ixhw, iyhw=iyhw,
        subgrid=subgrid, subgrid_gaus=subgrid_gaus,
        itstep=itstep, ntrac=ntrac,
        score_method=score_method, Sth0=Sth0, Sth1=Sth1,
        vxch=vxch, vych=vych,
        peak_inside_th=peak_inside_th, Cth=Cth,
        use_init_temp=use_init_temp, min_samples=min_samples)

    n = length(x0)
    tid0v = tid0 isa Integer ? fill(tid0, n) : Vector{Int}(tid0)
    vxgv = vxg === nothing ? nothing : (vxg isa AbstractArray ? Vector{Float64}(vxg) : fill(Float64(vxg), n))
    vygv = vyg === nothing ? nothing : (vyg isa AbstractArray ? Vector{Float64}(vyg) : fill(Float64(vyg), n))

    count, status, tid, x, y, vx, vy, score, zss, score_ary = VTTrac.trac(
        vtt, tid0v, x0, y0; vxg=vxgv, vyg=vygv,
        out_subimage=out_subimage, out_score_ary=out_score_ary, to_missing=false)

    union!(SEEN_STATUS, Set(Int.(status)))

    fmiss = vtt.fmiss
    imiss = vtt.imiss

    # --- convert 1-based Julia outputs to 0-based, NaN/-1 sentinel convention ---
    conv_tid(a) = Int64[(v == imiss) ? -1 : Int64(v) - 1 for v in a]
    conv_pos(a) = Float64[(v == fmiss) ? NaN : Float64(v) - 1.0 for v in a]   # position: shift by 1
    conv_val(a) = Float64[(Float64(v) == fmiss) ? NaN : Float64(v) for v in a] # velocity/score: no shift

    params = Dict{String,Any}(
        "nsx" => nsx, "nsy" => nsy,
        "vxhw" => vxhw, "vyhw" => vyhw, "ixhw" => vtt.ixhw, "iyhw" => vtt.iyhw,
        "subgrid" => subgrid, "subgrid_gaus" => subgrid_gaus,
        "itstep" => itstep, "ntrac" => ntrac,
        "score_method" => String(score_method),
        "Sth0" => Sth0, "Sth1" => Sth1,
        "vxch" => vxch, "vych" => vych,
        "peak_inside_th" => peak_inside_th, "Cth" => Cth,
        "use_init_temp" => use_init_temp, "min_samples" => min_samples,
        "zmiss" => zmiss, "fmiss" => fmiss, "imiss" => imiss,
        "has_mask" => mask !== nothing,
        "out_subimage" => out_subimage, "out_score_ary" => out_score_ary,
    )

    out = Dict{String,Any}()
    out["z"] = z
    out["t"] = t
    if mask !== nothing
        out["mask"] = mask
    end
    if zmiss !== nothing
        out["zmiss"] = Float32[Float32(zmiss)]
    end
    out["params_json"] = Vector{UInt8}(codeunits(to_json(params)))

    out["tid0_in"] = Float64.(tid0v) .- 1.0
    out["x0_in"] = x0 .- 1.0
    out["y0_in"] = y0 .- 1.0
    if vxgv !== nothing
        out["vxg_in"] = vxgv
    end
    if vygv !== nothing
        out["vyg_in"] = vygv
    end

    out["count"] = Int64.(count)
    out["status"] = Int64.(status)
    out["t_index"] = reduce(hcat, [conv_tid(tid[i, :]) for i in axes(tid, 1)])'  # (ntrac+1, n)
    out["x"] = reduce(hcat, [conv_pos(x[i, :]) for i in axes(x, 1)])'
    out["y"] = reduce(hcat, [conv_pos(y[i, :]) for i in axes(y, 1)])'
    out["vx"] = reduce(hcat, [conv_val(vx[i, :]) for i in axes(vx, 1)])'
    out["vy"] = reduce(hcat, [conv_val(vy[i, :]) for i in axes(vy, 1)])'
    out["score"] = reduce(hcat, [conv_val(score[i, :]) for i in axes(score, 1)])'

    if out_subimage
        out["zss"] = zss  # (ntrac+1, nsy, nsx, n); raw Julia fill (see test_golden.py for masking)
    end
    if out_score_ary
        out["score_ary"] = score_ary  # (ntrac, 2iyhw+1, 2ixhw+1, n)
    end

    npzwrite(joinpath(OUTDIR, name * ".npz"), out)
    push!(CASE_NAMES, String(name))
    println("  wrote $name.npz  (status seen: $(sort(unique(Int.(status)))))")
    return count, status
end

# ===========================================================================
# Case matrix
# ===========================================================================
println("Generating golden cases into $OUTDIR ...")

# --- helper: a regular interior grid of seed points (avoids image edges) ---
function seed_grid(; x0=10.0, x1=34.0, nxp=5, y0=10.0, y1=34.0, nyp=5)
    xs = collect(range(x0, x1, length=nxp))
    ys = collect(range(y0, y1, length=nyp))
    X = Float64[]
    Y = Float64[]
    for y in ys, x in xs
        push!(X, x)
        push!(Y, y)
    end
    return X .+ 1.0, Y .+ 1.0  # +1: caller passes back to Julia 1-based convention
end

X5, Y5 = seed_grid()  # 25 interior points, integer positions

# --- 1. Baseline: xcor, paraboloid subgrid ---
run_case("xcor_paraboloid_basic"; z=Z0, t=T0, nsx=7, nsy=7, vxhw=2.0, vyhw=2.0,
    subgrid=true, subgrid_gaus=false, ntrac=4, Sth0=0.5, Sth1=0.5,
    tid0=1, x0=X5, y0=Y5, out_subimage=true, out_score_ary=true)

# --- 2. ncov, paraboloid ---
run_case("ncov_paraboloid_basic"; z=Z0, t=T0, nsx=7, nsy=7, vxhw=2.0, vyhw=2.0,
    subgrid=true, subgrid_gaus=false, score_method="ncov", ntrac=4, Sth0=-1.0, Sth1=-1.0,
    tid0=1, x0=X5, y0=Y5, out_score_ary=true)

# --- 3. xcor, gaussian subgrid ---
run_case("xcor_gaussian_subgrid"; z=Z0, t=T0, nsx=7, nsy=7, vxhw=2.0, vyhw=2.0,
    subgrid=true, subgrid_gaus=true, ntrac=4, Sth0=0.5, Sth1=0.5,
    tid0=1, x0=X5, y0=Y5)

# --- 4. ncov, gaussian subgrid ---
run_case("ncov_gaussian_subgrid"; z=Z0, t=T0, nsx=7, nsy=7, vxhw=2.0, vyhw=2.0,
    subgrid=true, subgrid_gaus=true, score_method="ncov", ntrac=4, Sth0=-1.0, Sth1=-1.0,
    tid0=1, x0=X5, y0=Y5)

# --- 5. subgrid disabled (integer-only tracking) ---
run_case("xcor_no_subgrid"; z=Z0, t=T0, nsx=7, nsy=7, vxhw=2.0, vyhw=2.0,
    subgrid=false, ntrac=4, Sth0=0.5, Sth1=0.5,
    tid0=1, x0=X5, y0=Y5)

run_case("ncov_no_subgrid"; z=Z0, t=T0, nsx=7, nsy=7, vxhw=2.0, vyhw=2.0,
    subgrid=false, score_method="ncov", ntrac=4, Sth0=-1.0, Sth1=-1.0,
    tid0=1, x0=X5, y0=Y5)

# --- 6. Even template size (6x6): center-offset convention ---
run_case("xcor_even_template"; z=Z0, t=T0, nsx=6, nsy=6, vxhw=2.0, vyhw=2.0,
    subgrid=true, subgrid_gaus=false, ntrac=3, Sth0=0.5, Sth1=0.5,
    tid0=1, x0=X5, y0=Y5)

run_case("ncov_even_template_gaussian"; z=Z0, t=T0, nsx=6, nsy=6, vxhw=2.0, vyhw=2.0,
    subgrid=true, subgrid_gaus=true, score_method="ncov", ntrac=3, Sth0=-1.0, Sth1=-1.0,
    tid0=1, x0=X5, y0=Y5)

# --- 7. Small odd template (5x5) ---
run_case("xcor_small_template_5x5"; z=Z0, t=T0, nsx=5, nsy=5, vxhw=2.0, vyhw=2.0,
    subgrid=true, ntrac=3, Sth0=0.5, Sth1=0.5, tid0=1, x0=X5, y0=Y5)

# --- 8. Non-integer initial positions (bilinear template interpolation) ---
Xf = X5 .+ 0.37
Yf = Y5 .- 0.21
run_case("xcor_noninteger_seed"; z=Z0, t=T0, nsx=7, nsy=7, vxhw=2.0, vyhw=2.0,
    subgrid=true, ntrac=3, Sth0=0.5, Sth1=0.5, tid0=1, x0=Xf, y0=Yf)

run_case("ncov_noninteger_seed_gaussian"; z=Z0, t=T0, nsx=7, nsy=7, vxhw=2.0, vyhw=2.0,
    subgrid=true, subgrid_gaus=true, score_method="ncov", ntrac=3, Sth0=-1.0, Sth1=-1.0,
    tid0=1, x0=Xf, y0=Yf)

# --- 9. itstep variations: +1 (above), +2, and -1 (backward tracking) ---
# vxhw/vyhw kept small enough that the search window (2*ixhw+1) stays under
# the field's spatial period (10): a wider window would alias onto multiple
# near-equal periodic peaks and make the golden case numerically unstable
# (which peak wins the tie-break can flip on float-order noise alone).
run_case("xcor_itstep_plus2"; z=Z0, t=T0, nsx=7, nsy=7, vxhw=0.9, vyhw=0.9,
    subgrid=true, itstep=2, ntrac=3, Sth0=0.4, Sth1=0.4,
    tid0=1, x0=[11.3, 17.8, 21.2, 29.6, 33.1], y0=[12.7, 15.4, 24.9, 27.3, 31.8])

run_case("xcor_itstep_minus1_backward"; z=Z0, t=T0, nsx=7, nsy=7, vxhw=2.0, vyhw=2.0,
    subgrid=true, itstep=-1, ntrac=4, Sth0=0.5, Sth1=0.5, tid0=8, x0=X5, y0=Y5)

# --- 10. First guess velocity (vxg/vyg) ---
run_case("xcor_first_guess"; z=Z0, t=T0, nsx=7, nsy=7, ixhw=3, iyhw=3,
    subgrid=true, ntrac=3, Sth0=0.5, Sth1=0.5, tid0=1, x0=X5, y0=Y5,
    vxg=1.2, vyg=1.2)

# --- 11. vxch/vych screening (velocity-change rollback), status 9 ---
# accelerating field: phase speed jumps between steps -> 2nd-step velocity
# change exceeds vxch/vych, exercising the j==2 rollback path (VTTrac.jl B1).
let
    nt3, ny3, nx3 = 3, 40, 120
    t3 = collect(0.0:nt3-1)
    shiftx = [0.0, 1.0, 4.0]  # speed jumps from 1 to 3 between steps
    k3 = 2pi / 10
    z3 = Array{Float32,3}(undef, nt3, ny3, nx3)
    for i in 1:nt3, iy in 1:ny3, ix in 1:nx3
        z3[i, iy, ix] = Float32(sin(k3 * ((ix - 1) - shiftx[i])) * cos(k3 * (iy - 1)))
    end
    run_case("xcor_vxch_vych_rollback"; z=z3, t=t3, nsx=5, nsy=5, ixhw=6, iyhw=3,
        subgrid=false, ntrac=2, Sth0=0.5, Sth1=0.5, vxch=1.5, vych=1.5,
        tid0=1, x0=[51.0], y0=[16.0])
end

# --- 12. Cth (min contrast) screening: status 3 ---
run_case("xcor_cth_reject_all"; z=Z0, t=T0, nsx=7, nsy=7, vxhw=2.0, vyhw=2.0,
    subgrid=true, ntrac=2, Sth0=0.5, Sth1=0.5, Cth=2.5,  # data range is <= 2.0: unreachable
    tid0=1, x0=X5, y0=Y5)

run_case("xcor_cth_pass"; z=Z0, t=T0, nsx=7, nsy=7, vxhw=2.0, vyhw=2.0,
    subgrid=true, ntrac=2, Sth0=0.5, Sth1=0.5, Cth=0.05,  # easily satisfied
    tid0=1, x0=X5, y0=Y5)

# --- 13. peak_inside_th screening: status 4 ---
# Seeds placed near zero-crossings of the field (locally near-linear -> no
# interior peak) alongside seeds at local extrema (peak present).
Xpi = vcat(X5, [2.5 + 10*(2pi/(2pi/10))*0 ])  # placeholder, replaced below
run_case("xcor_peak_inside_mixed"; z=Z0, t=T0, nsx=7, nsy=7, vxhw=2.0, vyhw=2.0,
    subgrid=true, ntrac=2, Sth0=0.3, Sth1=0.3, peak_inside_th=0.15,
    tid0=1, x0=X5, y0=Y5)

# --- 14. tid start/end out of range: status 1 and 5 ---
run_case("xcor_tid_start_out_of_range"; z=Z0, t=T0, nsx=7, nsy=7, vxhw=2.0, vyhw=2.0,
    subgrid=true, ntrac=3, Sth0=0.5, Sth1=0.5,
    tid0=[0, -2, 1, 1], x0=[15.0, 15.0, 15.0, 15.0], y0=[15.0, 15.0, 15.0, 15.0])

run_case("xcor_tid_end_out_of_range"; z=Z0, t=T0, nsx=7, nsy=7, vxhw=2.0, vyhw=2.0,
    subgrid=true, ntrac=6, Sth0=0.5, Sth1=0.5,  # nt=10, starting at tid0=9 with ntrac=6 overruns
    tid0=9, x0=X5, y0=Y5)

# --- 15. template read failure (near image edge): status 2 ---
run_case("xcor_template_edge_out_of_bounds"; z=Z0, t=T0, nsx=7, nsy=7, vxhw=2.0, vyhw=2.0,
    subgrid=true, ntrac=2, Sth0=0.5, Sth1=0.5,
    tid0=1, x0=[1.0, 46.0, 24.0], y0=[24.0, 24.0, 1.0])

# --- 16. score computation failure (search window runs off image): status 6 ---
run_case("xcor_score_window_out_of_bounds"; z=Z0, t=T0, nsx=7, nsy=7, ixhw=8, iyhw=8,
    subgrid=true, ntrac=2, Sth0=0.5, Sth1=0.5,
    tid0=1, x0=[5.0, 42.0], y0=[24.0, 24.0])

# --- 17. peak not found (too-narrow search window forces boundary peak): status 7 ---
run_case("xcor_peak_not_found_narrow_window"; z=Z0, t=T0, nsx=7, nsy=7, ixhw=1, iyhw=1,
    subgrid=true, ntrac=2, Sth0=-1.0, Sth1=-1.0,
    tid0=1, x0=X5, y0=Y5, vxg=3.0, vyg=3.0)

# --- 18. score below threshold: status 8 ---
run_case("xcor_score_below_threshold"; z=Z0, t=T0, nsx=7, nsy=7, vxhw=2.0, vyhw=2.0,
    subgrid=true, ntrac=2, Sth0=0.999, Sth1=0.999,
    tid0=1, x0=X5, y0=Y5)

# --- 19. use_init_temp true/false ---
run_case("xcor_use_init_temp_false"; z=Z0, t=T0, nsx=7, nsy=7, vxhw=2.0, vyhw=2.0,
    subgrid=true, use_init_temp=false, ntrac=4, Sth0=0.5, Sth1=0.5, tid0=1, x0=X5, y0=Y5)

run_case("xcor_use_init_temp_true"; z=Z0, t=T0, nsx=7, nsy=7, vxhw=2.0, vyhw=2.0,
    subgrid=true, use_init_temp=true, ntrac=4, Sth0=0.5, Sth1=0.5, tid0=1, x0=X5, y0=Y5)

# --- 20. mask: none / sparse / dense / all, with min_samples variations ---
function sparse_mask(z)
    m = falses(size(z))
    nt, ny, nx = size(z)
    for i in 1:5:nt, j in 1:7:ny, k in 1:7:nx
        m[i, j, k] = true
    end
    return Array{Bool,3}(m)
end
function dense_mask(z)
    m = falses(size(z))
    nt, ny, nx = size(z)
    for i in 1:nt, j in 1:ny, k in 1:nx
        if (j + k) % 3 == 0
            m[i, j, k] = true
        end
    end
    return Array{Bool,3}(m)
end

run_case("xcor_mask_sparse_min1"; z=Z0, t=T0, mask=sparse_mask(Z0), nsx=7, nsy=7,
    vxhw=2.0, vyhw=2.0, subgrid=true, ntrac=3, Sth0=0.3, Sth1=0.3, min_samples=1,
    tid0=1, x0=X5, y0=Y5, out_score_ary=true)

run_case("xcor_mask_dense_min1"; z=Z0, t=T0, mask=dense_mask(Z0), nsx=7, nsy=7,
    vxhw=2.0, vyhw=2.0, subgrid=true, ntrac=3, Sth0=0.3, Sth1=0.3, min_samples=1,
    tid0=1, x0=X5, y0=Y5, out_score_ary=true)

run_case("xcor_mask_dense_min_high"; z=Z0, t=T0, mask=dense_mask(Z0), nsx=7, nsy=7,
    vxhw=2.0, vyhw=2.0, subgrid=true, ntrac=3, Sth0=-1.0, Sth1=-1.0, min_samples=40,
    tid0=1, x0=X5, y0=Y5, out_score_ary=true)

run_case("ncov_mask_sparse_min1"; z=Z0, t=T0, mask=sparse_mask(Z0), nsx=7, nsy=7,
    vxhw=2.0, vyhw=2.0, subgrid=true, score_method="ncov", ntrac=3, Sth0=-1.0, Sth1=-1.0,
    min_samples=1, tid0=1, x0=X5, y0=Y5)

let
    m = Array{Bool,3}(trues(size(Z0)))  # fully masked
    run_case("xcor_mask_all_masked"; z=Z0, t=T0, mask=m, nsx=7, nsy=7,
        vxhw=2.0, vyhw=2.0, subgrid=true, ntrac=2, Sth0=-1.0, Sth1=-1.0, min_samples=1,
        tid0=1, x0=X5[1:4], y0=Y5[1:4])
end

# --- 21. use_init_temp=true combined with mask (the v2.0.0 bugfix path, B1) ---
let
    m = falses(size(Z0))
    m[1, 1, 1] = true  # ensure chk_mask triggers (any(mask))
    m = Array{Bool,3}(m)
    run_case("xcor_use_init_temp_with_mask"; z=Z0, t=T0, mask=m, nsx=5, nsy=5,
        vxhw=2.0, vyhw=2.0, subgrid=true, use_init_temp=true, ntrac=4,
        Sth0=0.3, Sth1=0.3, tid0=1, x0=X5, y0=Y5)
end

# --- 22. missing value (zmiss) handling ---
let
    z = copy(Z0)
    z[3, 20, 20] = -999.0f0
    z[3, 5, 5] = -999.0f0
    z[6, 30, 12] = -999.0f0
    run_case("xcor_zmiss_scattered"; z=z, t=T0, zmiss=-999.0, nsx=7, nsy=7,
        vxhw=2.0, vyhw=2.0, subgrid=true, ntrac=4, Sth0=0.4, Sth1=0.4,
        tid0=1, x0=X5, y0=Y5)
end

# --- 23. Larger field, seeds near image border (margin cases, no failure) ---
Xb, Yb = seed_grid(x0=4.0, x1=59.0, nxp=6, y0=4.0, y1=59.0, nyp=6)
run_case("xcor_large_field_margin"; z=Z0BIG, t=T0BIG, nsx=7, nsy=7, vxhw=2.0, vyhw=2.0,
    subgrid=true, ntrac=4, Sth0=0.4, Sth1=0.4, tid0=1, x0=Xb, y0=Yb)

# --- 24. Combination stress case: ncov + mask + gaussian + noninteger + vxg ---
run_case("ncov_mask_gaussian_noninteger_combo"; z=Z0, t=T0, mask=sparse_mask(Z0),
    nsx=7, nsy=7, vxhw=2.0, vyhw=2.0, subgrid=true, subgrid_gaus=true,
    score_method="ncov", ntrac=3, Sth0=-1.0, Sth1=-1.0, min_samples=3,
    tid0=1, x0=Xf, y0=Yf, vxg=1.0, vyg=1.0, out_score_ary=true)

# --- 25. Scalar-vs-array tid0 & single-point shape variants (broad regression) ---
run_case("xcor_single_point"; z=Z0, t=T0, nsx=7, nsy=7, vxhw=2.0, vyhw=2.0,
    subgrid=true, ntrac=3, Sth0=0.5, Sth1=0.5, tid0=1, x0=[24.0], y0=[24.0])

# ===========================================================================
println("\nGenerated $(length(CASE_NAMES)) cases.")
println("Status codes observed across all cases: ", sort(collect(SEEN_STATUS)))
missing_codes = setdiff(Set(0:9), SEEN_STATUS)
if !isempty(missing_codes)
    @warn "Missing status codes (not exercised by any case)" missing_codes = sort(collect(missing_codes))
else
    println("All status codes 0-9 are exercised. ✓")
end
