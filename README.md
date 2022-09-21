# PLUMED Masterclass 21.3: Umbrella sampling

## Origin

This masterclass was authored by Giovanni Bussi on February 15, 2021.

## Aims

In this Masterclass, we will discuss how to perform and analyze umbrella sampling simulations. We will learn how to introduce a bias potential with PLUMED, how to compute the free energy landscapes, and how to reweight the resulting ensembles. We will also learn how to compute statistical errors on the computed quantities.

## Objectives

Once you have completed this Masterclass you will be able to:

- Use PLUMED to run simulations using static bias potentials with different functional forms.
- Use WHAM to combine multiple simulations performed with different bias potentials.
- Reweight the resulting ensembles so as to obtain the free-energy profile as a function of a different variable.
- Calculate error bars on free energies and populations.
 
## Setting up PLUMED

If you have not yet set up PLUMED, you can find information about installing it [here](https://www.plumed.org/doc-v2.7/user-doc/html/masterclass-21-1.html). 

Once you have installed PLUMED, you will need to install GROMACS as well. In particular, you will need a special version of GROMACS that has been patched with PLUMED. You can obtain it by using conda with the following command

````
conda install --strict-channel-priority -c plumed/label/masterclass -c conda-forge gromacs
````

The `--strict-channel-priority` is necessary as it prevents your conda install from downloading packages from the `bioconda` channel. `bioconda` contains a version of GROMACS that is **not** patched with PLUMED and would not work here.

On Linux, the command above should install the following packages:

````
  gromacs            plumed/label/masterclass/linux-64::gromacs-2019.6-h3fd9d12_0
  libclang           conda-forge/linux-64::libclang-11.0.1-default_ha53f305_1
  libevent           conda-forge/linux-64::libevent-2.1.10-hcdb4288_3
  libhwloc           conda-forge/linux-64::libhwloc-1.11.13-h3c4fd83_0
  libllvm11          conda-forge/linux-64::libllvm11-11.0.1-hf817b99_0
  libpq              conda-forge/linux-64::libpq-12.3-h255efa7_3
  [ etc ... ]
````

The exact versions might be different. Notice, however, that GROMACS comes from the `plumed/label/masterclass` channel, whereas the required libraries come from the `conda-forge` channel.  To be sure the installed GROMACS is patched with PLUMED, try the following shell command:

````
gmx mdrun -h 2> /dev/null | grep -q plumed && echo ok
````

It should print `ok`.

Please ensure that you have setup PLUMED and GROMACS on your machine before starting the exercises. Also notice that in order to obtain a good performance it is better to compile GROMACS from source on the machine where you are running your simulations. You can find out from the PLUMED documention how to patch GROMACS with PLUMED so as to be able to install it from source. For this tutorial, the conda precompiled binaries will be sufficient.

## Resources

The data needed to execute the exercises of this Masterclass can be found on [GitHub](https://github.com/plumed/masterclass-21-3).
You can clone this repository locally on your machine using the following command:

````
git clone https://github.com/plumed/masterclass-21-3.git
````

The files you need for the exercise are in the folder called `data`.  You will find the following files in that folder:

- `topolA.tpr`: an input file that can be used to run a GROMACS simulation of alanine dipeptide starting from one 
  of the two main free-energy minima
- `topolB.tpr`: same as `topolA.tpr`, but starting from the other minimum.
- `wham.py`: a python script that can be used to perform binless WHAM analysis
- `reference.pdb`: a pdb file of alanine dipeptide (the system under study) that can be used as input for your MOLINFO commands

Notice that PLUMED input files have not been provided in the GitHub repository.  You must prepare these input files yourself using the templates below.

We would recommend that you run each exercise in separate sub-directories inside the root directory `masterclass-21-3`.

_All the exercises were tested with PLUMED version 2.7.0 and GROMACS 2019.6_

## Exercises

Throughout this tutorial we will run simulations of alanine dipeptide in vacuum using GROMACS and PLUMED.  Whereas this system is too simple to be considered a proper benchmark for enhanced sampling methods, it is complex enough to be used when learning about them. Some of the commands below are specific for GROMACS, but all the PLUMED input files are compatible with other MD engines as well.

### Exercise 1: Running an unbiased simulation with PLUMED

In other masterclasses you learned how to analyze trajectories a posteriori. One of the nice features of PLUMED is that the very same analysis can be done on the fly. In other words, you might compute your collective variables while GROMACS is running instead of waiting for the simulation to end. This can could be convenient if you want to run a large number of simulations, you already know what to compute, and you do not want to use too much disk space.

To run a simulation with GROMACS you have to type this command in the shell:

````
gmx mdrun -plumed plumed.dat -s topolA.tpr -nsteps 200000 -x traj_unbiased.xtc
````

Notice that the file `topolA.tpr` contains all the relevant information (simulation parameters, initial conditions, etc.).  In this tutorial we will just need to play with the number of steps (200000 in the example above) and we will tune the name of the trajectory saved by GROMACS (traj_unbiased.xtc in the example above).

Also consider that GROMACS and PLUMED will never delete your files.  They will take backups when you try to overwrite them. GROMCAS backup files start with `#`, whereas PLUMED backup files start with `bck.`. Both GROMACS and PLUMED will complain if you try to create too many backups of the same file. It is thus recommended to regularly clean your directory with a command such as `rm -f \#* bck.*`

The command above will only succeed if a file named `plumed.dat` exists in the current directory.  You should already know to create such file from having completed other (more basic) tutorials. Those other tutorials should also have taught you how to compute histograms. You should be thus able to complete the template below and put it in a file named `plumed.dat`.

```plumed
#SOLUTIONFILE=work/plumed_ex1.dat
# vim:ft=plumed
MOLINFO STRUCTURE=reference.pdb
phi: TORSION ATOMS=__FILL__ # use MOLINFO shortcuts to identify the phi angle of the second residue
psi: TORSION ATOMS=__FILL__ # use MOLINFO shortcuts to identify the psi angle of the second residue

# use the command below to compute the histogram of phi
# we use a smooth kernel to produce a nicer graph here
# notice that when asking for numbers PLUMED is happy to accept strings such as "pi" meaning 3.14...
# also arithmetics is allowed (e.g. 2*pi could be used if necessary)
hhphi: HISTOGRAM ARG=__FILL__ STRIDE=100 GRID_MIN=-pi GRID_MAX=__FILL__ GRID_BIN=600 BANDWIDTH=0.05
ffphi: CONVERT_TO_FES GRID=hhphi # no need to set TEMP here, PLUMED will obtain it from GROMACS
DUMPGRID GRID=__FILL__ FILE=fes_phi.dat STRIDE=200000 # stride is needed here since PLUMED does not know when the simulation is over

# now add three more lines to compute and dump the free energy as a function of **psi** on a file names fes_psi.dat
__FILL__




PRINT __FILL__ # use this command to write phi and psi on a file named colvar.dat, every 100 steps
```

You can then monitor what happened during the simulation using the python script below.

```python

import plumed
import matplotlib.pyplot as plt
import numpy as np

# plot the time series of phi and psi
colvar=plumed.read_as_pandas("colvar.dat")
plt.plot(colvar.time,colvar.phi,"x",label="phi")
plt.plot(colvar.time,colvar.psi,"x",label="psi")
plt.xlabel("time")
plt.ylabel("$\phi$")
plt.legend()
plt.show()

# scatter plot with phi and psi
plt.plot(colvar.phi,colvar.psi,"x")
plt.xlabel("$\phi$")
plt.ylabel("$\psi$")
plt.xlim((-np.pi,np.pi))
plt.ylim((-np.pi,np.pi))
plt.show()

# FES as a function of phi
# we remove infinite and nans here
fes_phi=plumed.read_as_pandas("fes_phi").replace([np.inf, -np.inf], np.nan).dropna()
plt.plot(fes_phi.phi,fes_phi.ffphi)
plt.xlim((-np.pi,np.pi))
plt.xlabel("$\phi$")
plt.ylabel("$F(\phi)$")
plt.show()

# FES as a function of psi
# we remove infinite and nans here
fes_psi=plumed.read_as_pandas("fes_psi").replace([np.inf, -np.inf], np.nan).dropna()
plt.plot(fes_psi.psi,fes_psi.ffpsi)
plt.xlim((-np.pi,np.pi))
plt.xlabel("$\psi$")
plt.ylabel("$F(\psi)$")
plt.show()
```

In this manner you will be able to see (a) the time series of $\phi$ and $\psi$, (b) which portion of the $(\phi,\psi)$ space has been explored (c) the free energy as a function of $\phi$ and (d) the free energy as a function of $\psi$.  In particular, you should be able to identify two minima separated by a moderate barrier. Several transitions between these minima will be visible in the simulated time scale.

### Exercise 2: Running biased simulations with PLUMED

We are now ready to use PLUMED to perform the task it was originally designed for: biasing a simulation on the fly.  We will first try to add simple bias potentials that change the balance between the two minima we have obtained in the previous exercise.  Biasing is not particularly useful here since both minima were sampled anyways. However, it is instructive since we will be able to test the types
of analysis we did in [masterclass-21-2](https://plumed-school.github.io/lessons/21/002/data/). The two minima observed in the previous exercise are located at $\phi$ approximately equal to -2.5 and -1.5and are separated by a barrier located at $\phi$ approximately equal to -2.  We can predict that adding a bias potential in the form `-A*sin(x+2)`, with $A$ being positive and large enough, should favor the minimum at -1.5 at the expense of the other.

We can try this type of calculation with A=10 by using an input file like the following one:

```plumed
#SOLUTIONFILE=work/plumed_ex2.dat
# vim:ft=plumed
__FILL__ # compute phi and psi here, as in the previous input file

# fill in with the required function (i.e. -10*sin(phi+2)):
f: CUSTOM ARG=phi FUNC=__FILL__ PERIODIC=NO 

# this command allows to add a bias potential equal to f
BIASVALUE ARG=f

# here you can paste the same HISTOGRAM/CONVERT_TO_FES/DUMPGRID commands that you used in
# the previous exercise. let's just write the results on different files
hhphi: __FILL__
ffphi: __FILL__
DUMPGRID FILE=fes_phi_biased1.dat __FILL__
hhpsi: __FILL__
ffpsi: __FILL__
DUMPGRID FILE=fes_psi_biased1.dat __FILL__

lw: REWEIGHT_BIAS

# and here we do the same again, but this time using LOGWEIGHTS.
# these free energies will be printed on files fes_phi_biased1r.dat and 
# fes_psi_biased1r.dat and will be reweighted so as to be unbiased
hhphir: HISTOGRAM __FILL__ LOGWEIGHTS=lw
ffphir: __FILL__
DUMPGRID FILE=fes_phi_biased1r.dat __FILL__

hhpsir: HISTOGRAM __FILL__ LOGWEIGHTS=lw
ffpsir: __FILL__
DUMPGRID FILE=fes_psi_biased1r.dat __FILL__

PRINT __FILL__ # monitor what's happening, as before, writing on file plumed_colvar1.dat
```

Call this file `plumed_biased1.dat` and run the simulation using this command:

````
gmx mdrun -plumed plumed_biased1.dat -s topolA.tpr -nsteps 200000 -x traj_comp_biased1.xtc
````

Notice that we are storing the trajectory on a separate file. We will need both the unbiased and the biased trajectory later.

You can then monitor what happened during the simulation using these commands

```python
colvar=plumed.read_as_pandas("colvar_biased1.dat")
plt.plot(colvar.time,colvar.phi,"x",label="phi")
plt.plot(colvar.time,colvar.psi,"x",label="psi")
plt.xlabel("time")
plt.ylabel("$\phi$")
plt.legend()
plt.show()

plt.plot(colvar.phi,colvar.psi,"x")
plt.xlabel("$\phi$")
plt.ylabel("$\psi$")
plt.xlim((-np.pi,np.pi))
plt.ylim((-np.pi,np.pi))
plt.show()

fes_phi=plumed.read_as_pandas("fes_phi.dat").replace([np.inf, -np.inf], np.nan).dropna()
plt.plot(fes_phi.phi,fes_phi.ffphi,label="original")
fes_phib=plumed.read_as_pandas("fes_phi_biased1.dat").replace([np.inf, -np.inf], np.nan).dropna()
plt.plot(fes_phib.phi,fes_phib.ffphi,label="biased")
fes_phir=plumed.read_as_pandas("fes_phi_biased1r.dat").replace([np.inf, -np.inf], np.nan).dropna()
plt.plot(fes_phir.phi,fes_phir.ffphir,label="reweighted")
plt.legend()
plt.xlim((-np.pi,np.pi))
plt.xlabel("$\phi$")
plt.ylabel("$F(\phi)$")
plt.show()

fes_psi=plumed.read_as_pandas("fes_psi.dat").replace([np.inf, -np.inf], np.nan).dropna()
plt.plot(fes_psi.psi,fes_psi.ffpsi,label="original")
fes_psib=plumed.read_as_pandas("fes_psi_biased1.dat").replace([np.inf, -np.inf], np.nan).dropna()
plt.plot(fes_psib.psi,fes_psib.ffpsi,label="biased")
fes_psir=plumed.read_as_pandas("fes_psi_biased1r.dat").replace([np.inf, -np.inf], np.nan).dropna()
plt.plot(fes_psir.psi,fes_psir.ffpsir,label="reweighted")
plt.legend()
plt.xlim((-np.pi,np.pi))
plt.xlabel("$\psi$")
plt.ylabel("$F(\psi)$")
plt.show()
```

The free-energy plots here will include three lines: (a) the original one (as obtained from the previous exercise), (b) the one sampled in this second simulation, where the minimum at higher phi appears more stable, and (c) the reweighted free-energy from the second simulation. This final free energy surface should look closer to the original one.  Notice that the relationship between plots (b) and (c) is straightforward when the analyzed variable is phi: the third line could have been obtained by just adding $-10 \sin(x+2)$  to the second line.  However, when analyzing psi the effect is less obvious.

### Exercise 3: Combining statistics from biased and unbiased simulations

We will now run a new simulation that is even more biased by choosing the prefactor A=20.  This will further favour the minimum at higher values of phi.  Please prepare a `plumed_biased2.dat` input file that is identical to `plumed_biased1.dat` except that:

- it applies a bias with -20*sin(phi+2)
- all the written files are named with a `2` suffix (e.g. `fes_phi_biased2.dat`).

Then run the new simulation using the following command:

````
gmx mdrun -plumed plumed_biased2.dat -s topolA.tpr -nsteps 200000 -x traj_comp_biased2.xtc
````

You can analyze this new simulation as you did the previous one.  You should notice that the peak at higher $\phi$ is even more favored now.

We will now use WHAM to combine the three simulations we have done so far, namely:

1. unbiased (`traj_comp_unbiased.xtc`),
2. biased (`traj_comp_biased1.xtc`), and
3. more biased (`traj_comp_biased2.xtc`).

The first thing we have to do is to concatenate the three files usinig:

````
gmx trjcat -cat -f traj_comp_unbiased.xtc traj_comp_biased1.xtc traj_comp_biased2.xtc -o traj_comp_cat.xtc
````

Notice that this new trajectory contains all the simulated frames, and we __do not__ explicitly track which of the simulation each frame originated from. 

We then have to compute the bias that would have been felt in each of the three runs above by each of the frames in the concatenated trajectory. We could do this reusing the three plumed input files we used above, but it is perhaps clearer to do it with a single input file in this case. You can use an input like this one:

```plumed
#SOLUTIONFILE=work/plumed_ex3.dat
# vim:ft=plumed
MOLINFO STRUCTURE=reference.pdb
phi: TORSION __FILL__ 
psi: TORSION __FILL__

# here we list the bias potentials that we used in the three simulations we want to combine
b0: CUSTOM ARG=phi FUNC=0.0 PERIODIC=NO
b1: CUSTOM ARG=phi FUNC=__FILL__ PERIODIC=NO # fill here with the bias you used in the first biased simulation
b2: CUSTOM ARG=phi FUNC=__FILL__ PERIODIC=NO # fill here with the bias you used in the second biased simulation
PRINT ARG=phi,psi,b0,b1,b2 FILE=biases.dat
```

If you call your input file `plumed_wham.dat` you can then use the following command

````
plumed driver --plumed plumed_wham.dat --ixtc traj_comp_cat.xtc
````

We are now ready to use binless WHAM to compute the weights associated to the concatenated trajectory.  This can be done using the following script

```python
import wham
# this is to check that you actually imported the module located in this directory
# it should print /current/directory/wham.py
print(wham.__file__)

bias=plumed.read_as_pandas("biases.dat")
kBT=300*8.314462618*0.001 # use kJ/mol here
w=wham.wham(np.stack((bias.b0,bias.b1,bias.b2)).T,T=kBT)
print(w)

# you can then add a new column to the bias dataframe:
print(bias)
bias["logweights"]=w["logW"]
print(bias)

# and write it on disk
plumed.write_pandas(bias,"bias_wham.dat")
```

The procedure above computes the weights that can be associated to the concatenated trajectory so as to recover an unbiased distribution. Notice that, although we did not explicitly enforce it,
the weights only depend on the value of $\phi$ as only $\phi$ was biased in our simulations.  Thus, if we plot the weight as a function of $\phi$ all the points will collapse on a line.
This will not happen if we plot the weights as a function of $\psi$:

```python
plt.plot(bias.phi,bias.logweights,"x")
plt.show()
plt.plot(bias.psi,bias.logweights,"x")
plt.show()
```

Also notice the form of the weights as a function of $\phi$. This graph is telling us that points at larger $\phi$ are weighted less. This makes sense, since in two simulations out of three we were biasing the ensemble so as to visit those points **more** frequently than in the Boltzmann ensemble. However, the curve is not exactly a sin function. This precise curve can only be obtained by solving the WHAM equations self-consistently.

We can then read the obtained weights and use them in an analysis similar to the ones we have done above

```plumed
#SOLUTIONFILE=work/plumed_ex4.dat
# vim:ft=plumed
phi: READ FILE=bias_wham.dat VALUES=__FILL__
psi: READ FILE=bias_wham.dat VALUES=__FILL__
lw: READ FILE=bias_wham.dat VALUES=__FILL__

hhphi: HISTOGRAM ARG=phi GRID_MIN=-pi GRID_MAX=pi GRID_BIN=600 BANDWIDTH=0.05
ffphi: CONVERT_TO_FES GRID=hhphi
DUMPGRID GRID=ffphi FILE=fes_phi_cat.dat

hhpsi: HISTOGRAM ARG=psi GRID_MIN=-pi GRID_MAX=pi GRID_BIN=600 BANDWIDTH=0.05
ffpsi: CONVERT_TO_FES GRID=hhpsi
DUMPGRID GRID=ffpsi FILE=fes_psi_cat.dat

# we use a smooth kernel to produce a nicer graph here
hhphir: HISTOGRAM ARG=phi GRID_MIN=-pi GRID_MAX=pi GRID_BIN=600 BANDWIDTH=0.05 LOGWEIGHTS=lw
ffphir: CONVERT_TO_FES GRID=hhphir
DUMPGRID GRID=ffphir FILE=fes_phi_catr.dat

hhpsir: HISTOGRAM ARG=psi GRID_MIN=-pi GRID_MAX=pi GRID_BIN=600 BANDWIDTH=0.05 LOGWEIGHTS=lw
ffpsir: CONVERT_TO_FES GRID=hhpsir
DUMPGRID GRID=ffpsir FILE=fes_psi_catr.dat
```

Use the template above to produce a file names `plumed_wham.dat` and run the following command:

````
plumed driver --noatoms --plumed plumed_wham.dat --kt 2.4943387854
````

You can now show the result by using the following script

```python
colvar=plumed.read_as_pandas("bias_wham.dat")
plt.plot(colvar.time,colvar.phi,"x",label="phi")
plt.plot(colvar.time,colvar.psi,"x",label="psi")
plt.xlabel("time")
plt.ylabel("$\phi$")
plt.legend()
plt.show()

plt.plot(colvar.phi,colvar.psi,"x")
plt.xlabel("$\phi$")
plt.ylabel("$\psi$")
plt.xlim((-np.pi,np.pi))
plt.ylim((-np.pi,np.pi))
plt.show()

fes_phi=plumed.read_as_pandas("fes_phi.dat").replace([np.inf, -np.inf], np.nan).dropna()
plt.plot(fes_phi.phi,fes_phi.ffphi,label="original")
fes_phib=plumed.read_as_pandas("fes_phi_cat.dat").replace([np.inf, -np.inf], np.nan).dropna()
plt.plot(fes_phib.phi,fes_phib.ffphi,label="biased")
fes_phir=plumed.read_as_pandas("fes_phi_catr.dat").replace([np.inf, -np.inf], np.nan).dropna()
plt.plot(fes_phir.phi,fes_phir.ffphir,label="reweighted")
plt.legend()
plt.xlim((-np.pi,np.pi))
plt.xlabel("$\phi$")
plt.ylabel("$F(\phi)$")
plt.show()

fes_psi=plumed.read_as_pandas("fes_psi.dat").replace([np.inf, -np.inf], np.nan).dropna()
plt.plot(fes_psi.psi,fes_psi.ffpsi,label="original")
fes_psib=plumed.read_as_pandas("fes_psi_cat.dat").replace([np.inf, -np.inf], np.nan).dropna()
plt.plot(fes_psib.psi,fes_psib.ffpsi,label="biased")
fes_psir=plumed.read_as_pandas("fes_psi_catr.dat").replace([np.inf, -np.inf], np.nan).dropna()
plt.plot(fes_psir.psi,fes_psir.ffpsir,label="reweighted")
plt.legend()
plt.xlim((-np.pi,np.pi))
plt.xlabel("$\psi$")
plt.ylabel("$F(\psi)$")
plt.show()
```

When you look at the first graph (i.e., the concatenated time series) you should notice that as you proceed towards the end the frames are coming from the most biased simulation, and thus will have a systematically larger value of $\phi$. The reweighted ensemble will be quite similar to the original one in the end.

### Exercise 4: Enhancing conformational transitions with multiple-windows umbrella sampling

So far we only considered the fast transition between the two minima at negative values of $\phi$. However, alanine dipeptide has another metastable state at positive $\phi$, that is separated by the minima we have seen by a large barrier.  In line of principle, one could use the trick above to discount the effect of any arbitrary biasing potential.  For instance, if one had an estimate of a relevant free-energy barrier, a bias potential approximately equal to the negative of the barrier would make the residual free energy basically barrierless. However, this procedure requires knowing the free energy profile in advance. As we will see in the next Masterclass, metadynamics provides exactly a self-healing method that adjusts the bias until dynamics becomes diffusive.

In the context of the static potentials that we are using here, this is a very difficult approach.  In the tutorial [lugano-2](https://www.plumed.org/doc-v2.7/user-doc/html/lugano-2.html) you can find an attempt to derive a potential that cancels the large barrier between the states.  When this bias is used the number of transitions between the important metastable states increases dramatically. In this class we will take a more pragmatic approach and directly introduce the method that is commonly used in these situations, namely "multiple windows umbrella sampling".

To do so we will have to run a number of separate simulations, each of them with a bias potential designed to maintain the $\phi$ variable in a narrow range. If the consecutive windows overlap sufficiently, and if the windows span the relevant region in the $\phi$ variable, we will be able to reconstruct the full free-energy profile.

We will use 32 simulations, with `RESTRAINT` potentials centered at uniformly spaced values of $\phi$ in the range between -pi and pi.  In order to avoid the repetition of the first potential, recalling that $\phi$ is periodic, we can generate the list with the command

```python
import numpy as np
at=np.linspace(-np.pi,np.pi,32,endpoint=False)
print(at)
```

Now you should write a python script that generate 32 plumed input files (you can call them `plumed_0.dat`, `plumed_1.dat`, etc)
that look like this one

```plumed
#SOLUTIONFILE=work/plumed_ex5.dat
# vim:ft=plumed
MOLINFO STRUCTURE=reference.pdb
phi: TORSION ATOMS=__FILL__
psi: TORSION ATOMS=__FILL__
bb: RESTRAINT ARG=phi KAPPA=200.0 AT=__FILL__
PRINT ARG=phi,psi,bb.bias FILE=__FILE__ STRIDE=100
# make sure that each simulation writes on a different file
# e.g. colvar_multi_0.dat, colvar_multi_1.dat, ...
```

Now run the 32 simulations with commands like these ones

````
gmx mdrun -plumed plumed_0.dat -s topolA.tpr -nsteps 200000 -x traj_comp_0.xtc
gmx mdrun -plumed plumed_1.dat -s topolA.tpr -nsteps 200000 -x traj_comp_1.xtc
# etc.
````

As we did before, you should now concatenate the resulting trajectories

````
gmx trjcat -cat -f traj_comp_[0-9].xtc traj_comp_[0-9][0-9].xtc -o traj_multi_cat.xtc
````

and analyze them with plumed driver:

````
plumed driver --plumed plumed_0.dat --ixtc traj_multi_cat.xtc --trajectory-stride 100
plumed driver --plumed plumed_1.dat --ixtc traj_multi_cat.xtc --trajectory-stride 100
# etc.
````

Notice the `--trajectory-stride` option: we are telling the driver that the trajectory was saved every 100 frames.  In this manner, we will be able to recycle the same plumed input file that we used while running the simulations.  This approach is more robust than the one we used above, where a single plumed file was used to reproduce all the biases, and is thus recommended in production cases.
Also notice that it is common practice to ignore the initial part of each simulation, since it typically contains a transient that would make the final result less robust. For simplicitly, we are not doing it here.  Finally, notice that in the last exercise we wrote a new input file just for writing the three biases. Here instea we recycled the input files used to run the 32 simulations.

We can now read the produced files and process them with our WHAM script

```python
col=[]
for i in range(32):
    col.append(plumed.read_as_pandas("colvar_multi_" + str(i)+".dat"))
# notice that this is the concatenation of 32 trajectories with 2001 frames each
    plt.plot(col[i].phi[2001*i:2001*(i+1)],col[i].psi[2001*i:2001*(i+1)],"x")
plt.xlabel("$\phi$")
plt.ylabel("$\psi$")
plt.show()
# in this graph you can appreciate which region was sampled by each simulation

bias=np.zeros((len(col[0]["bb.bias"]),32))
for i in range(32):
    bias[:,i]=col[i]["bb.bias"][-len(bias):]
w=wham.wham(bias,T=kBT)
plt.plot(w["logW"])
plt.show()
colvar=col[0]
colvar["logweights"]=w["logW"]
plumed.write_pandas(colvar,"bias_multi.dat")
```

Now you can do something similar to what we did before, namely:

- write a file to be read by PLUMED.
- process it with plumed driver so as to compute the free energy with respect to $\phi$ and $\psi$
- plot the free energy as a function of $\phi$ and $\psi$

How do the resulting free-energy profiles compares with the one we obtained before?

### Exercise 5: Computing populations and errors

Now that you have a simulation that sampled all the relevant metastable states, you can analyze them in many different ways.  In particular, you can compute the average of any quantity just using the weights that we obtained.  This can be done directly with PLUMED, as we have seen in other masterclasses, by using the [AVERAGE](https://www.plumed.org/doc-master/user-doc/html/_a_v_e_r_a_g_e.html) command.
However, we will use Python to compute averages here as this will make it easier to compute their statistical errors with bootstrap.

To complete this exercise you should compute the average value and the statistical error of the following quantities:

- The population of the metastable state B defined as $0<\phi<2$. Notice that the exact definition is not very important if the minimum is clearly delimited by high free-energy barriers (i.e., low probability regions).
- Free-energy difference between the state B defined as $0<\phi<2$ and the state A defined as $-pi<\phi<0$.  Notice that the free-energy difference is defined as $-k_B T \log\left( \frac{P_B}{P_A}\right)$, where $P_B$ is the population of state B and $P_A$ the population of state A.
- Average value of $\phi$ in state A and in state B. Warning: be careful computing the average of an angle?  Remember this quantity is periodic and that must be taken into account.

For each of these quantities, you should also compute the standard error. It is here recommended to compute the error using boostrap. You can do so by adjusting the following script

```python
# number of blocks
NB=10
# we reshape the bias so that it appears as NB times frames of each traj times number of biases
bb=bias.reshape((32,-1,32))[:,-2000:,:].reshape((32,NB,2000//NB,32))
# we reshape the trajectory so that it appears as NB times frames of each traj times number of biases
cc=np.array(col[0].phi).reshape((32,-1))[:,-2000:].reshape((32,NB,2000//NB))

# we first analyse the complete trajectory:
tr=cc.flatten()
is_in_B=np.int_(np.logical_and(tr>0,tr<2))
w0=wham.wham(bb.reshape((-1,32)),T=kBT)

print("population:",np.average(is_in_B,weights=np.exp(w0["logW"])))
pop=[]
for i in range(200):
    # we then analyze the bootstrapped trajectories
    print(i)
    c=np.random.choice(NB,NB)
    w=wham.wham(bb[:,c,:,:].reshape((-1,32)),T=kBT)
    tr=cc[:,c,:].flatten()
    is_in_B=np.int_(np.logical_and(tr>0,tr<2))
    pop.append(np.average(is_in_B,weights=np.exp(w["logW"])))
```

### Exercise 6: Effect of initial conditions

Repeat the second to last exercise above using `topolB.tpr` as an input file. This will make the simulations start from the minimum in B.  Plot the free energy as a function of $\phi$ and as a function of $\psi$. How does the result you obtain differ from the one you obtained when you performed the earlier exercise?
