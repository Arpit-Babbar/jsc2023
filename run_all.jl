using Pkg
using Trixi: trixi_include
Pkg.add("Git")
run(`$(git()) clone https://github.com/arpit-babbar/Tenkai.jl`)
include("dmr/elixir_double_mach.jl")
include("dmr/run_double_mach_reflection.jl")
include("kh/elixir_euler_kelvin_helmholtz_instability.jl")
include("kh/run_kh_schaal.jl")
include("m2000/elixir_euler_astrophysical_m2000.jl")
include("m2000/run_jet_m2000.jl")
include("rp2d/elixir_euler_rp2.jl")
include("rp2d/run_rp2d_12.jl")
include("sedov/elixir_euler_sedov_blast_wave.jl")
include("sedov/run_sedov_blast2d.jl")
