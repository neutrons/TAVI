using Sunny, GLMakie, LinearAlgebra

units = Units(:meV, :angstrom);

a = 5.45;
latvecs = lattice_vectors(a, a, a, 55.13, 55.13, 55.13);
positions = [[0.35300, 0.35300, 0.35300]];
cryst = Crystal(latvecs, positions, 148, choice="R"; types=["Ni"]);
view_crystal(cryst)

print_symmetry_table(cryst, 8.0)

# published parameters
# Jab1 = -0.15;
# Jab2 = -0.05;
# Jc1 = 0;
# Jc2 = 0.27;
# Jab3 = -0.05;
# Jc3 = 0.27;
# D = 0.1; # (meV) hard-axis
#
Jab1 = -0.12;
Jc1 = 0.24;
Jab2 = -0.08;
Jc2 = 0.25;
Jab3 = -0.07;
Jc3 = 0.23;
D = 0.12; # (meV) hard-axis

sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole, seed=1)
set_exchange!(sys, Jab1, Bond(1, 2, [-1, 0, 0]));
set_exchange!(sys, Jab2, Bond(1, 1, [-1, 1, 0]));
set_exchange!(sys, Jc1, Bond(1, 2, [0, 0, 0]));
set_exchange!(sys, Jc2, Bond(1, 1, [1, 0, 0]));
set_exchange!(sys, Jab3, Bond(1, 2, [-1, -1, 1]));
set_exchange!(sys, Jc3, Bond(2, 1, [1, 1, 0]));
n = normalize(cryst.latvecs * [1, 1, 1]);
set_onsite_coupling!(sys, S -> D * (n' * S)^2, 1);


# find the ground state
sys = resize_supercell(sys, (2, 2, 2))
randomize_spins!(sys)
minimize_energy!(sys)
plot_spins(sys; color=[S[1] for S in sys.dipoles])
print_wrapped_intensities(sys)

suggest_magnetic_supercell([[1 / 2, 1 / 2, 1 / 2]])
sys_min = reshape_supercell(sys, [1 1 1; -1 1 0; 0 0 1])
randomize_spins!(sys_min)
minimize_energy!(sys_min);
plot_spins(sys_min; color=[S[1] for S in sys_min.dipoles], ghost_radius=12)

# Linear spin wave theory
formfactors = [1 => FormFactor("Ni2")]
measure = ssf_perp(sys_min; formfactors)
swt = SpinWaveTheory(sys_min; measure)

kernel = gaussian(fwhm=0.1)

qs = [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1], [1.5, 1.5, 1.5], [2, 2, 2], [2.5, 2.5, 2.5],];
path = q_space_path(cryst, qs, 2000)
energies = range(0, 4, 1000)
res = intensities(swt, path; energies, kernel)
plot_intensities(res; units, colormap=:jet, colorrange=(0, 20))


qs = [[0, 0, 0], [0.5, 0.5, 0.5], [1, 1, 1], [1.5, 1.5, 1.5], [2, 2, 2]]
path = q_space_path(cryst, qs, 200)
res = intensities_bands(swt, path)
plot_intensities(res; units, colorrange=(0, 5))


# grid = q_space_grid(cryst, [1, 1, 0], range(-1, 1, 200), [0, 0, 1], (0, 10); orthogonalize=true)
# res = intensities(swt, grid; energies=[3.7], kernel)
# plot_intensities(res; units, colormap=:jet, colorrange=(0, 20))
