#!/bin/bash

# parameters
box_size=150       # simulation box size (L)
num_of_atoms=1000  # number of atoms/beads (M)
ratio=2.0          # feedback ratio in Sneppen model (F = alpha / (1 - alpha))
bond_energy=1.0    # bond energy in pairwise potential (epsilon)
cut_off=2.5        # pairwise potential cut-off parameter (r_c)
frac_static=0.0    # fraction of bookmark/static atoms (phi)
cluster_size=0     # cluster size (n_c)
tcolour=10         # recolour time (in Brownian time units)
tmax=1000000       # maximum simulation time (in Brownian time units), must be multiple of tcolour 
teq=10000          # total equlibration steps
run=1              # trial number
atom_type=0        # configuration for non-static atoms: 0 = random, 1 = active, 2 = unmark, 3 = inactive, 4 = random uniform
static_type=random # configuration for static atoms (random, cluster, mixed)
dumpxyz=dump       # dump or nodump
dumpstate=state    # nostate or state
rundir="./"        # main directory where the code will be run
print_freq=1000    # dump atoms' positions frequency (in Brownian time units)

type_of_atoms=6  # types of atoms
delta_t=0.01     # time step size in Brownian time units
colour_step=$(bc <<< "$tcolour/$delta_t")
max_iter=$(bc <<< "$tmax/$delta_t/$colour_step")
print_freq=$(bc <<< "$print_freq/$delta_t")  # actual print frequency in simulation steps

# generate a random seed for initialising dna and running lammps
seed=$(python GetRandom.py 100000)

# make execution directory
rundir="${rundir}/run_${run}_f_${ratio}_e_${bond_energy}_phi_${frac_static}_nc_${cluster_size}"
mkdir $rundir

# copy execution code to run directory
cp -r dna_epigenetics $rundir

# generate the data/output file names
name="L_${box_size}_N_${num_of_atoms}_f_${ratio}_e_${bond_energy}_rc_${cut_off}_phi_${frac_static}_nc_${cluster_size}_t_${tmax}_run_${run}"
thermo_file="thermo_${name}.dat"
xyz_file="vmd_${name}.xyz"
init_file="init_${name}.in"
in_file="dna_${name}.in"
out_file="dna_${name}.out"
stats_file="stats_${name}.dat"
state_file="none"  
pos_file="none"

if [ $dumpstate != "nostate" ]; then
    state_file="state_${name}.dat"
fi

if [ $dumpxyz != "nodump" ]; then
    pos_file="pos_${name}.dat"
fi

# create the lammps command file based on template
lammps_file="epi_${name}.lam"
file="${rundir}/${lammps_file}"

# choose template depending on initial conditions (collasped/swollen)
cp epigenetics-swollen.lam $file

# replace macros in template with input values
sed -i -- "s/NUMOFATOMS/${num_of_atoms}/g" $file
sed -i -- "s/RATIO/${ratio}/g" $file
sed -i -- "s/BONDENERGY/${bond_energy}/g" $file
sed -i -- "s/CUTOFF/${cut_off}/g" $file
sed -i -- "s/MAXITER/${max_iter}/g" $file
sed -i -- "s/COLOURSTEP/${colour_step}/g" $file
sed -i -- "s/SEED/${seed}/g" $file
sed -i -- "s/RUN/${run}/g" $file
sed -i -- "s/XYZFILE/${xyz_file}/g" $file
sed -i -- "s/INITFILE/${init_file}/g" $file
sed -i -- "s/INFILE/${in_file}/g" $file
sed -i -- "s/OUTFILE/${out_file}/g" $file
sed -i -- "s/LAMMPSFILE/${lammps_file}/g" $file
sed -i -- "s/STATEFILE/${state_file}/g" $file
sed -i -- "s/STATSFILE/${stats_file}/g" $file
sed -i -- "s/POSFILE/${pos_file}/g" $file
sed -i -- "s/PRINTFREQ/${print_freq}/g" $file
sed -i -- "s/DELTAT/${delta_t}/g" $file

# calculate equilibrium times
t_soft=10 # equilibrate time with soft potential (in Brownian time units)
equil_soft=$(bc <<< "$t_soft/$delta_t")
equil_total=$(bc <<< "$teq/$delta_t")
equil_fene=$(bc <<< "$equil_total-$equil_soft")

sed -i -- "s/EQUIL_SOFT/${equil_soft}/g" $file
sed -i -- "s/EQUIL_FENE/${equil_fene}/g" $file
sed -i -- "s/EQUIL_TOTAL/${equil_total}/g" $file


# initialise dna strand (ordered/disordered)
seed2=$(bc <<< "$seed+34987")

java dna_epigenetics.LAMMPSIO $num_of_atoms $type_of_atoms $box_size $box_size $box_size $seed2 $atom_type $static_type $frac_static $cluster_size "${rundir}/${init_file}"

# clear any previous entries in the state, stats, and pos file
if [ $state_file != "none" ]; then
    > "${rundir}/${state_file}"
fi

if [ $stats_file != "none" ]; then
    > "${rundir}/${stats_file}"
fi

if [ $pos_file != "none" ]; then
    > "${rundir}/${pos_file}"
fi
