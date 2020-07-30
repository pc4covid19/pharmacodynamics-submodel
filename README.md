# COVID19 pharmacodynamics-submodel
This model simulates the COVID19 pharmacodynamics response. See https://github.com/pc4covid19/COVID19 for full model information.\
To model the single-cell pharmacodynamics, we first modeled the rate of drug internalization as diffusion through the cell membrane, allowing us to track the accumulated total drug in each cell. Then, we model the total effect E for each cell based on the internal drug concentration using a Hill function. The effect is then used to modulate a vector r of cell parameters from their untreated rates r_0 to their rates at maximum drug efficacy r_max when E  = Emax. In this investigation, r includes the ACE2 dynamics, viral replication and export rate. In addition, we can model the impact of partial drug efficacy by multiplying the effect by a sensitivity S, where S = 0 denotes a drug with no efficacy (e.g., it fails to bind its target), and S = 100% is an idea drug that fully binds and inhibits its target. 

You can try to run our cloud-hosted nanoHUB app for fun! \
https://nanohub.org/tools/virion2pbpd

* * * 

## Release summary:
### 0.1.0:
This is initial release. 


