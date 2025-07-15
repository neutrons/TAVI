from juliacall import Main as jl

jl.include("./scripts/convolution/SW01_FM_Heseinberg_chain.jl")


print(jl.cryst)

# qs = [[0, 0, 0], [1, 0, 0]]
# path = jl.q_space_path(jl.cryst, qs, 400)

# res = jl.intensities_bands(jl.swt, path)
