using Sunny, LinearAlgebra

function model()
    a = 5.45
    latvecs = lattice_vectors(a, a, a, 55.13, 55.13, 55.13)
    positions = [[0.35300, 0.35300, 0.35300]]
    cryst = Crystal(latvecs, positions, 148, choice="R"; types=["Ni"])

    sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole, seed=7)
    # published parameters
    # Jab1 = -0.15;
    # Jab2 = -0.05;
    # Jc1 = 0;
    # Jc2 = 0.27;
    # Jab3 = -0.05;
    # Jc3 = 0.27;
    # D = 0.1; # (meV) hard-axis
    # 
    Jab1 = -0.12
    Jc1 = 0.24
    Jab2 = -0.08
    Jc2 = 0.25
    Jab3 = -0.07
    Jc3 = 0.23
    D = 0.12 # (meV) hard-axis

    sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole, seed=7)
    set_exchange!(sys, Jab1, Bond(1, 2, [-1, 0, 0]))
    set_exchange!(sys, Jab2, Bond(1, 1, [-1, 1, 0]))
    set_exchange!(sys, Jc1, Bond(1, 2, [0, 0, 0]))
    set_exchange!(sys, Jc2, Bond(1, 1, [1, 0, 0]))
    set_exchange!(sys, Jab3, Bond(1, 2, [-1, -1, 1]))
    set_exchange!(sys, Jc3, Bond(2, 1, [1, 1, 0]))
    n = normalize(cryst.latvecs * [1, 1, 1])
    set_onsite_coupling!(sys, S -> D * (n' * S)^2, 1)

    sys_min = reshape_supercell(sys, [1 1 1; -1 1 0; 0 0 1])
    randomize_spins!(sys_min)
    minimize_energy!(sys_min)

    # Linear spin wave theory
    formfactors = [1 => FormFactor("Ni2")]
    measure = ssf_perp(sys_min; formfactors)
    swt = SpinWaveTheory(sys_min; measure)
    return swt
end


# function disp(swt, qpts)
#     res = dispersion(swt, qpts)
#     return res
# end

function disp(swt, qpts)
    res = intensities_bands(swt, qpts)
    return res.disp
end

function inten(swt, qpts)
    res = intensities_bands(swt, qpts)
    return res.data
end

# units = Units(:meV, :angstrom);
# swt = model()
#qs = [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1]]
#path = q_space_path(swt.sys.crystal, qs, 200)
# res = intensities_bands(swt, path)
# using GLMakie
# plot_intensities(res; units, colorrange=(0, 5))