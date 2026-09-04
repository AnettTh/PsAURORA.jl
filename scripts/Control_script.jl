using AURORA

## Setting parameters
altitude_lims = [100, 600];     # (km) altitude limits of the ionosphere
θ_lims = 180:-30:0              # (°) angle-limits for the electron beams
E_max = 15000;                  # (eV) upper limit to the energy grid
B_angle_to_zenith = 13;         # (°) angle between the B-field line and the zenith

msis_file = find_msis_file(
    year=2005, month=10, day=8, hour=22, minute=0, lat=70, lon=19, height=85:1:700
    );
iri_file = find_iri_file(
    year=2005, month=10, day=8, hour=22, minute=0, lat=70, lon=19, height=85:1:700
    );

## Build the model
model = AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file, B_angle_to_zenith)

## Define where to save the results
output = AuroraOutputManager("data/boris_test"; overwrite=false)

## Define input flux
flux = InputFlux(FlatSpectrum(1e-2; E_min=14500), SinusoidalFlickering(5.0);
                 beams=1, z_source=4*RE / 1e3, propagation=:fieldline)


## Create and run the simulation
mode = TimeDependent(duration = 5,           # (s) total simulation time
                     dt = 0.005,                # (s) time step for saving data
                     CFL_number = 128,
                     #n_loop = 2,             # (optional) define manually the number of loops to run
                     #max_memory_gb = 2.0,     # (optional) or determine n_loop based on limit memory usage
                     )

sim = AuroraSimulation(model, flux, output; mode)

##
run!(sim)

## Run the analysis
make_Ie_top_file(sim)
make_volume_excitation_file(sim)
make_current_file(sim)
make_column_excitation_file(sim)
