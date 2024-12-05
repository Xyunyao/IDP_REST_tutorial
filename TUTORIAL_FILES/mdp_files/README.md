# mdp parameter files used for GROMACS simulation

The files provided here will be used in various steps as an input for `gmx grompp`. 
    - `minimz.mdp` : Used to generate `.tpr` for energy minimization and can also be used for ion generation step.
    - `NVT.mdp` : Used to generate `.tpr` for thermally heating up the solvated protein box and generate velocities. We use velocity rescale thermostat as a temperature bath.
    - `NPT0.mdp` : Used to generate `.tpr` for quick equilibration to get pressure around 1 bar. We use Berendsen barostat.
    - `NPT1.mdp` : Used to generate `.tpr` to perform optimal pressure equilibration. We can use Parrinello-Rahman/C-rescale barostat.
    - `prod.mdp` : Used to generate `.tpr` to perform REST2 production run. 