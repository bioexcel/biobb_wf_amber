# Protein MD Setup tutorial using BioExcel Building Blocks (biobb) 
### --***AmberTools package version***--
**Based on the [MDWeb](http://mmb.irbbarcelona.org/MDWeb2/) [Amber FULL MD Setup tutorial](https://mmb.irbbarcelona.org/MDWeb2/help.php?id=workflows#AmberWorkflowFULL)**

***
This tutorial aims to illustrate the process of **setting up a simulation system** containing a **protein**, step by step, using the **BioExcel Building Blocks library (biobb)** wrapping the **AmberTools** utility from the **AMBER package**. The particular example used is the **Lysozyme** protein (PDB code 1AKI, [https://doi.org/10.2210/pdb1AKI/pdb](https://doi.org/10.2210/pdb1AKI/pdb)). 
***

## Settings

### Biobb modules used

 - [biobb_io](https://github.com/bioexcel/biobb_io): Tools to fetch biomolecular data from public databases.
 - [biobb_amber](https://github.com/bioexcel/biobb_amber): Tools to setup and run Molecular Dynamics simulations with AmberTools.
 - [biobb_analysis](https://github.com/bioexcel/biobb_analysis): Tools to analyse Molecular Dynamics trajectories.
 
### Auxiliary libraries used

* [jupyter](https://jupyter.org/): Free software, open standards, and web services for interactive computing across all programming languages.
* [plotly](https://plot.ly/python/offline/): Python interactive graphing library integrated in Jupyter notebooks.
* [nglview](https://nglviewer.org/#nglview): Jupyter/IPython widget to interactively view molecular structures and trajectories in notebooks.
* [simpletraj](https://github.com/arose/simpletraj): Lightweight coordinate-only trajectory reader based on code from GROMACS, MDAnalysis and VMD.
* [gfortran](https://anaconda.org/conda-forge/gfortran): Fortran 95/2003/2008/2018 compiler for GCC, the GNU Compiler Collection.
* [libgfortran5](https://anaconda.org/conda-forge/libgfortran5): Fortran compiler and libraries from the GNU Compiler Collection.

### Conda Installation

```console
git clone https://github.com/bioexcel/biobb_wf_amber_md_setup.git
cd biobb_wf_amber_md_setup
conda env create -f conda_env/environment.yml
conda activate biobb_wf_amber_md_setup
jupyter-notebook biobb_wf_amber_md_setup/notebooks/md_setup/biobb_wf_amber_md_setup.ipynb
```

***
## Pipeline steps
 1. [Input Parameters](#input)
 2. [Fetching PDB Structure](#fetch)
 3. [Preparing PDB file for AMBER](#pdb4amber)
 4. [Create Protein System Topology](#top)
 5. [Energetically Minimize the Structure](#minv)
 6. [Create Solvent Box and Solvating the System](#box)
 7. [Adding Ions](#ions)
 8. [Energetically Minimize the System](#min)
 9. [Heating the System](#heating)
 10. [Equilibrate the System (NVT)](#nvt)
 11. [Equilibrate the System (NPT)](#npt)
 12. [Free Molecular Dynamics Simulation](#free)
 13. [Post-processing and Visualizing Resulting 3D Trajectory](#post)
 14. [Output Files](#output)
 15. [Questions & Comments](#questions)
 
***
<img src="https://bioexcel.eu/wp-content/uploads/2019/04/Bioexcell_logo_1080px_transp.png" alt="Bioexcel2 logo"
	title="Bioexcel2 logo" width="400" />
***


## Initializing colab
The two cells below are used only in case this notebook is executed via **Google Colab**. Take into account that, for running conda on **Google Colab**, the **condacolab** library must be installed. As [explained here](https://pypi.org/project/condacolab/), the installation requires a **kernel restart**, so when running this notebook in **Google Colab**, don't run all cells until this **installation** is properly **finished** and the **kernel** has **restarted**.


```python
# Only executed when using google colab
import sys
if 'google.colab' in sys.modules:
  import subprocess
  from pathlib import Path
  try:
    subprocess.run(["conda", "-V"], check=True)
  except FileNotFoundError:
    subprocess.run([sys.executable, "-m", "pip", "install", "condacolab"], check=True)
    import condacolab
    condacolab.install()
    # Clone repository
    repo_URL = "https://github.com/bioexcel/biobb_wf_amber_md_setup.git"
    repo_name = Path(repo_URL).name.split('.')[0]
    if not Path(repo_name).exists():
      subprocess.run(["mamba", "install", "-y", "git"], check=True)
      subprocess.run(["git", "clone", repo_URL], check=True)
      print("⏬ Repository properly cloned.")
    # Install environment
    print("⏳ Creating environment...")
    env_file_path = f"{repo_name}/conda_env/environment.yml"
    subprocess.run(["mamba", "env", "update", "-n", "base", "-f", env_file_path], check=True)
    print("🎨 Install NGLView dependencies...")
    subprocess.run(["mamba", "install", "-y", "-c", "conda-forge", "nglview==3.0.8", "ipywidgets=7.7.2"], check=True)
    print("👍 Conda environment successfully created and updated.")
```


```python
# Enable widgets for colab
if 'google.colab' in sys.modules:
  from google.colab import output
  output.enable_custom_widget_manager()
  # Change working dir
  import os
  os.chdir("biobb_wf_amber_md_setup/biobb_wf_amber_md_setup/notebooks/md_setup")
  print(f"📂 New working directory: {os.getcwd()}")
```

<a id="input"></a>
## Input parameters
**Input parameters** needed:
 - **pdbCode**: PDB code of the protein structure (e.g. 1AKI, [https://doi.org/10.2210/pdb1AKI/pdb](https://doi.org/10.2210/pdb1AKI/pdb))


```python
import nglview
import ipywidgets
import plotly
from plotly import subplots
import plotly.graph_objs as go

pdbCode = "1aki"
```

<a id="fetch"></a>
***
## Fetching PDB structure
Downloading **PDB structure** with the **protein molecule** from the RCSB PDB database.<br>
Alternatively, a **PDB file** can be used as starting structure. <br>
***
**Building Blocks** used:
 - [pdb](https://biobb-io.readthedocs.io/en/latest/api.html#module-api.pdb) from **biobb_io.api.pdb**
***


```python
# Import module
from biobb_io.api.pdb import pdb

# Create properties dict and inputs/outputs
downloaded_pdb = pdbCode+'.pdb'

prop = {
    'pdb_code': pdbCode
}

#Create and launch bb
pdb(output_pdb_path=downloaded_pdb,
    properties=prop)
```

<a id="vis3D"></a>
### Visualizing 3D structure
Visualizing the downloaded/given **PDB structure** using **NGL**:    


```python
# Show protein
view = nglview.show_structure_file(downloaded_pdb)
view.add_representation(repr_type='ball+stick', selection='all')
view._remote_call('setSize', target='Widget', args=['','600px'])
view
```

<img src='_static/md_setup/ngl1.png'></img>

<a id="pdb4amber"></a>
***
## Preparing PDB file for AMBER
Before starting a **protein MD setup**, it is always strongly recommended to take a look at the initial structure and try to identify important **properties** and also possible **issues**. These properties and issues can be serious, as for example the definition of **disulfide bridges**, the presence of a **non-standard aminoacids** or **ligands**, or **missing residues**. Other **properties** and **issues** might not be so serious, but they still need to be addressed before starting the **MD setup process**. **Missing hydrogen atoms**, presence of **alternate atomic location indicators** or **inserted residue codes** (see [PDB file format specification](https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM)) are examples of these not so crucial characteristics. Please visit the [AMBER tutorial: Building Protein Systems in Explicit Solvent](http://ambermd.org/tutorials/basic/tutorial7/index.php) for more examples. **AmberTools** utilities from **AMBER MD package** contain a tool able to analyse **PDB files** and clean them for further usage, especially with the **AmberTools LEaP program**: the **pdb4amber tool**. The next step of the workflow is running this tool to analyse our **input PDB structure**.<br>

For the particular **Lysosyme** example, the most important property that is identified by the **pdb4amber** utility is the presence of **disulfide bridges** in the structure. Those are marked changing the residue names **from CYS to CYX**, which is the code that **AMBER force fields** use to distinguish between cysteines forming or not forming **disulfide bridges**. This will be used in the following step to correctly form a **bond** between these cysteine residues. 

We invite you to check what the tool does with different, more complex structures (e.g. PDB code [6N3V](https://www.rcsb.org/structure/6N3V)). 

***
**Building Blocks** used:
 - [pdb4amber_run](https://biobb-amber.readthedocs.io/en/latest/pdb4amber.html#pdb4amber-pdb4amber-run-module) from **biobb_amber.pdb4amber.pdb4amber_run**
***


```python
# Import module
from biobb_amber.pdb4amber.pdb4amber_run import pdb4amber_run

# Create prop dict and inputs/outputs
output_pdb4amber_path = 'structure.pdb4amber.pdb'

# Create and launch bb
pdb4amber_run(  input_pdb_path=downloaded_pdb,
            output_pdb_path=output_pdb4amber_path)
```

<a id="top"></a>
***
## Create protein system topology
**Building AMBER topology** corresponding to the protein structure.<br>

*IMPORTANT: the previous pdb4amber building block is changing the proper cysteines residue naming in the PDB file from CYS to CYX so that this step can automatically identify and add the disulfide bonds to the system topology.*<br>

The **force field** used in this tutorial is [**ff14SB**](https://doi.org/10.1021/acs.jctc.5b00255), an evolution of the **ff99SB** force field with improved accuracy of protein side chains and backbone parameters. **Water** molecules type used in this tutorial is [**tip3p**](https://doi.org/10.1021/jp003020w).<br>
Adding **side chain atoms** and **hydrogen atoms** if missing. Forming **disulfide bridges** according to the info added in the previous step. <br>

Generating three output files: 
- **AMBER structure** (PDB file)
- **AMBER topology** (AMBER [Parmtop](https://ambermd.org/FileFormats.php#topology) file)
- **AMBER coordinates** (AMBER [Coordinate/Restart](https://ambermd.org/FileFormats.php#restart) file) 
***
**Building Blocks** used:
 - [leap_gen_top](https://biobb-amber.readthedocs.io/en/latest/leap.html#module-leap.leap_gen_top) from **biobb_amber.leap.leap_gen_top**
***


```python
# Import module
from biobb_amber.leap.leap_gen_top import leap_gen_top

# Create prop dict and inputs/outputs
output_pdb_path = 'structure.leap.pdb'
output_top_path = 'structure.leap.top'
output_crd_path = 'structure.leap.crd'

prop = {
    "forcefield" : ["protein.ff14SB"]
}

# Create and launch bb
leap_gen_top(input_pdb_path=output_pdb4amber_path,
           #input_pdb_path=downloaded_pdb,
           output_pdb_path=output_pdb_path,
           output_top_path=output_top_path,
           output_crd_path=output_crd_path,
           properties=prop)
```

### Visualizing 3D structure
Visualizing the **PDB structure** using **NGL**. <br>
Try to identify the differences between the structure generated for the **system topology** and the **original one** (e.g. hydrogen atoms).   


```python
# Show protein
view = nglview.show_structure_file(output_pdb_path)
view.add_representation(repr_type='ball+stick', selection='all')
view._remote_call('setSize', target='Widget', args=['','600px'])
view
```

<img src='_static/md_setup/ngl2.png'></img>

<a id="minv"></a>
## Energetically minimize the structure

**Energetically minimize** the **protein structure** (in vacuo) using the **sander tool** from the **AMBER MD package**. This step is **relaxing the structure**, usually **constrained**, especially when coming from an X-ray **crystal structure**. <br/>

The **miminization process** is done in two steps:
- [Step 1](#minv_1): **Hydrogen** minimization, applying **position restraints** (50 Kcal/mol.$Å^{2}$) to the **protein heavy atoms**.
- [Step 2](#minv_2): **System** minimization, applying **position restraints** (50 Kcal/mol.$Å^{2}$) to the **protein heavy atoms**.
***
**Building Blocks** used:
 - [sander_mdrun](https://biobb-amber.readthedocs.io/en/latest/sander.html#module-sander.sander_mdrun) from **biobb_amber.sander.sander_mdrun**
 - [process_minout](https://biobb-amber.readthedocs.io/en/latest/process.html#module-process.process_minout) from **biobb_amber.process.process_minout**
***

<a id="minv_1"></a>
### Step 1: Minimize Hydrogens
**Hydrogen** minimization, applying **position restraints** (50 Kcal/mol.$Å^{2}$) to the **protein heavy atoms**.


```python
# Import module
from biobb_amber.sander.sander_mdrun import sander_mdrun

# Create prop dict and inputs/outputs
output_h_min_traj_path = 'sander.h_min.x'
output_h_min_rst_path = 'sander.h_min.rst'
output_h_min_log_path = 'sander.h_min.log'

prop = {
    'simulation_type' : "min_vacuo",
    "mdin" : { 
        'maxcyc' : 500,
        'ntpr' : 5,
        'ntr' : 1,
        'restraintmask' : '\":*&!@H=\"',
        'restraint_wt' : 50.0
    }
}

# Create and launch bb
sander_mdrun(input_top_path=output_top_path,
            input_crd_path=output_crd_path,
            input_ref_path=output_crd_path,
            output_traj_path=output_h_min_traj_path,
            output_rst_path=output_h_min_rst_path,
            output_log_path=output_h_min_log_path,
            properties=prop)
```

### Checking Energy Minimization results
Checking **energy minimization** results. Plotting **potential energy** along time during the **minimization process**.


```python
# Import module
from biobb_amber.process.process_minout import process_minout

# Create prop dict and inputs/outputs
output_h_min_dat_path = 'sander.h_min.energy.dat'

prop = {
    "terms" : ['ENERGY']
}

# Create and launch bb
process_minout(input_log_path=output_h_min_log_path,
            output_dat_path=output_h_min_dat_path,
            properties=prop)  
```


```python
import plotly.graph_objs as go

with open(output_h_min_dat_path, 'r') as energy_file:
    x, y = zip(*[
        (float(line.split()[0]), float(line.split()[1]))
        for line in energy_file
        if not line.startswith(("#", "@"))
        if float(line.split()[1]) < 1000
    ])

# Create a scatter plot
fig = go.Figure(data=go.Scatter(x=x, y=y, mode='lines'))

# Update layout
fig.update_layout(title="Energy Minimization",
                  xaxis_title="Energy Minimization Step",
                  yaxis_title="Potential Energy kcal/mol",
                  height=600)

# Show the figure (renderer changes for colab and jupyter)
rend = 'colab' if 'google.colab' in sys.modules else ''
fig.show(renderer=rend)
```

<img src='_static/md_setup/plot01.png'></img>

<a id="minv_2"></a>
### Step 2: Minimize the system
**System** minimization, applying **position restraints** (50 Kcal/mol.$Å^{2}$) to the **protein heavy atoms**.


```python
# Import module
from biobb_amber.sander.sander_mdrun import sander_mdrun

# Create prop dict and inputs/outputs
output_n_min_traj_path = 'sander.n_min.x'
output_n_min_rst_path = 'sander.n_min.rst'
output_n_min_log_path = 'sander.n_min.log'

prop = {
    'simulation_type' : "min_vacuo",
    "mdin" : { 
        'maxcyc' : 500,
        'ntpr' : 5,
        'ntr' : 1,
        'restraintmask' : '\":*&!@H=\"',
        'restraint_wt' : 50.0
    }
}

# Create and launch bb
sander_mdrun(input_top_path=output_top_path,
            input_crd_path=output_h_min_rst_path,
            input_ref_path=output_h_min_rst_path,
            output_traj_path=output_n_min_traj_path,
            output_rst_path=output_n_min_rst_path,
            output_log_path=output_n_min_log_path,
            properties=prop)
```

### Checking Energy Minimization results
Checking **energy minimization** results. Plotting **potential energy** by time during the **minimization process**.


```python
# Import module
from biobb_amber.process.process_minout import process_minout

# Create prop dict and inputs/outputs
output_n_min_dat_path = 'sander.n_min.energy.dat'

prop = {
    "terms" : ['ENERGY']
}

# Create and launch bb
process_minout(input_log_path=output_n_min_log_path,
            output_dat_path=output_n_min_dat_path,
            properties=prop)
```


```python
import plotly.graph_objs as go

with open(output_n_min_dat_path, 'r') as energy_file:
    x, y = zip(*[
        (float(line.split()[0]), float(line.split()[1]))
        for line in energy_file
        if not line.startswith(("#", "@"))
        if float(line.split()[1]) < 1000
    ])

# Create a scatter plot
fig = go.Figure(data=go.Scatter(x=x, y=y, mode='lines'))

# Update layout
fig.update_layout(title="Energy Minimization",
                  xaxis_title="Energy Minimization Step",
                  yaxis_title="Potential Energy kcal/mol",
                  height=600)

# Show the figure (renderer changes for colab and jupyter)
rend = 'colab' if 'google.colab' in sys.modules else ''
fig.show(renderer=rend)
```

<img src='_static/md_setup/plot02.png'></img>

<a id="box"></a>
***
## Create solvent box and solvating the system
Define the unit cell for the **protein structure MD system** to fill it with water molecules.<br>
A **truncated octahedron box** is used to define the unit cell, with a **distance from the protein to the box edge of 9.0 Angstroms**. <br>
The solvent type used is the default **TIP3P** water model, a generic 3-point solvent model.
***
**Building Blocks** used:
 - [amber_to_pdb](https://biobb-amber.readthedocs.io/en/latest/ambpdb.html#module-ambpdb.amber_to_pdb) from **biobb_amber.ambpdb.amber_to_pdb**
 - [leap_solvate](https://biobb-amber.readthedocs.io/en/latest/leap.html#module-leap.leap_solvate) from **biobb_amber.leap.leap_solvate** 
***

### Getting minimized structure
Getting the result of the **energetic minimization** and converting it to **PDB format** to be then used as input for the **water box generation**. <br/>This is achieved by converting from **AMBER topology + coordinates** files to a **PDB file** using the **ambpdb** tool from the **AMBER MD package**.


```python
# Import module
from biobb_amber.ambpdb.amber_to_pdb import amber_to_pdb

# Create prop dict and inputs/outputs
output_ambpdb_path = 'structure.ambpdb.pdb'

# Create and launch bb
amber_to_pdb(input_top_path=output_top_path,
            input_crd_path=output_n_min_rst_path,
            output_pdb_path=output_ambpdb_path
            )
```

### Create water box
Define the **unit cell** for the **protein structure MD system** and fill it with **water molecules**.<br/>


```python
# Import module
from biobb_amber.leap.leap_solvate import leap_solvate

# Create prop dict and inputs/outputs
output_solv_pdb_path = 'structure.solv.pdb'
output_solv_top_path = 'structure.solv.parmtop'
output_solv_crd_path = 'structure.solv.crd'

prop = {
    "forcefield" : ["protein.ff14SB"],
    "water_type": "TIP3PBOX",
    "distance_to_molecule": "9.0",   
    "box_type": "truncated_octahedron"
}

# Create and launch bb
leap_solvate(input_pdb_path=output_ambpdb_path,
           output_pdb_path=output_solv_pdb_path,
           output_top_path=output_solv_top_path,
           output_crd_path=output_solv_crd_path,
           properties=prop)
```

<a id="ions"></a>
## Adding ions

**Neutralizing** the system and adding an additional **ionic concentration** using the **leap tool** from the **AMBER MD package**. <br/>
Using **Sodium (Na+)** and **Chloride (Cl-)** counterions and an **additional ionic concentration** of 150mM.
***
**Building Blocks** used:
 - [leap_add_ions](https://biobb-amber.readthedocs.io/en/latest/leap.html#module-leap.leap_add_ions) from **biobb_amber.leap.leap_add_ions**
***


```python
# Import module
from biobb_amber.leap.leap_add_ions import leap_add_ions

# Create prop dict and inputs/outputs
output_ions_pdb_path = 'structure.ions.pdb'
output_ions_top_path = 'structure.ions.parmtop'
output_ions_crd_path = 'structure.ions.crd'

prop = {
    "forcefield" : ["protein.ff14SB"],
    "neutralise" : True,
    "positive_ions_type": "Na+",
    "negative_ions_type": "Cl-",
    "ionic_concentration" : 150, # 150mM
    "box_type": "truncated_octahedron"
}

# Create and launch bb
leap_add_ions(input_pdb_path=output_solv_pdb_path,
           output_pdb_path=output_ions_pdb_path,
           output_top_path=output_ions_top_path,
           output_crd_path=output_ions_crd_path,
           properties=prop)
```

### Visualizing 3D structure
Visualizing the **protein system** with the newly added **solvent box** and **counterions** using **NGL**.<br> Note the **truncated octahedron box** filled with **water molecules** surrounding the **protein structure**, as well as the randomly placed **positive** (Na+, blue) and **negative** (Cl-, gray) **counterions**.


```python
# Show protein
view = nglview.show_structure_file(output_ions_pdb_path)
view.clear_representations()
view.add_representation(repr_type='cartoon', selection='protein')
view.add_representation(repr_type='ball+stick', selection='solvent')
view.add_representation(repr_type='spacefill', selection='Cl- Na+')
view._remote_call('setSize', target='Widget', args=['','600px'])
view
```

<img src='_static/md_setup/ngl3.png'></img>

<a id="min"></a>
## Energetically minimize the system

**Energetically minimize** the **system** (protein structure + solvent + ions) using the **sander tool** from the **AMBER MD package**. **Restraining heavy atoms** with a force constant of 15 Kcal/mol.$Å^{2}$ to their initial positions.

- [Step 1](#emStep1): Energetically minimize the **system** through 500 minimization cycles.
- [Step 2](#emStep2): Checking **energy minimization** results. Plotting energy by time during the **minimization** process. 
***
**Building Blocks** used:
 - [sander_mdrun](https://biobb-amber.readthedocs.io/en/latest/sander.html#module-sander.sander_mdrun) from **biobb_amber.sander.sander_mdrun**
 - [process_minout](https://biobb-amber.readthedocs.io/en/latest/process.html#module-process.process_minout) from **biobb_amber.process.process_minout**
***

<a id="emStep1"></a>
### Step 1: Running Energy Minimization
The **minimization** type of the **simulation_type property** contains the main default parameters to run an **energy minimization**:

-  imin  = 1 ;    Minimization flag, perform an energy minimization.
-  maxcyc = 500;  The maximum number of cycles of minimization.
-  ntb = 1;       Periodic boundaries: constant volume.
-  ntmin = 2;     Minimization method: steepest descent.


In this particular example, the method used to run the **energy minimization** is the default **steepest descent**, with a **maximum number of 500 cycles** and **periodic conditions**. 


```python
# Import module
from biobb_amber.sander.sander_mdrun import sander_mdrun

# Create prop dict and inputs/outputs
output_min_traj_path = 'sander.min.x'
output_min_rst_path = 'sander.min.rst'
output_min_log_path = 'sander.min.log'

prop = {
    "simulation_type" : "minimization",
    "mdin" : { 
        'maxcyc' : 300, # Reducing the number of minimization steps for the sake of time
        'ntr' : 1,      # Overwritting restrain parameter
        'restraintmask' : '\"!:WAT,Cl-,Na+\"',      # Restraining solute
        'restraint_wt' : 50.0                       # With a force constant of 50 Kcal/mol*A2
    }
}

# Create and launch bb
sander_mdrun(input_top_path=output_ions_top_path,
            input_crd_path=output_ions_crd_path,
            input_ref_path=output_ions_crd_path,
            output_traj_path=output_min_traj_path,
            output_rst_path=output_min_rst_path,
            output_log_path=output_min_log_path,
            properties=prop)
```

<a id="emStep2"></a>
### Step 2: Checking Energy Minimization results
Checking **energy minimization** results. Plotting **potential energy** along time during the **minimization process**.


```python
# Import module
from biobb_amber.process.process_minout import process_minout

# Create prop dict and inputs/outputs
output_dat_path = 'sander.min.energy.dat'

prop = {
    "terms" : ['ENERGY']
}

# Create and launch bb
process_minout(input_log_path=output_min_log_path,
            output_dat_path=output_dat_path,
            properties=prop)
```


```python
import plotly.graph_objs as go

with open(output_dat_path, 'r') as energy_file:
    x, y = zip(*[
        (float(line.split()[0]), float(line.split()[1]))
        for line in energy_file
        if not line.startswith(("#", "@"))
        if float(line.split()[1]) < 1000
    ])

# Create a scatter plot
fig = go.Figure(data=go.Scatter(x=x, y=y, mode='lines'))

# Update layout
fig.update_layout(title="Energy Minimization",
                  xaxis_title="Energy Minimization Step",
                  yaxis_title="Potential Energy kcal/mol",
                  height=600)

# Show the figure (renderer changes for colab and jupyter)
rend = 'colab' if 'google.colab' in sys.modules else ''
fig.show(renderer=rend)
```

<img src='_static/md_setup/plot03.png'></img>

<a id="heating"></a>
## Heating the system

**Warming up** the **prepared system** using the **sander tool** from the **AMBER MD package**. Going from 0 to the desired **temperature**, in this particular example, 300K. **Solute atoms restrained** (force constant of 10 Kcal/mol). Length 5ps.
***
- [Step 1](#heatStep1): Warming up the **system** through 500 MD steps.
- [Step 2](#heatStep2): Checking results for the **system warming up**. Plotting **temperature** along time during the **heating** process. 
***
**Building Blocks** used:
 - [sander_mdrun](https://biobb-amber.readthedocs.io/en/latest/sander.html#module-sander.sander_mdrun) from **biobb_amber.sander.sander_mdrun**
 - [process_mdout](https://biobb-amber.readthedocs.io/en/latest/process.html#module-process.process_mdout) from **biobb_amber.process.process_mdout**
***

<a id="heatStep1"></a>
### Step 1: Warming up the system
The **heat** type of the **simulation_type property** contains the main default parameters to run a **system warming up**:

-  imin = 0;&nbsp;&nbsp;&nbsp;      Run MD (no minimization)
-  ntx = 5;&nbsp;&nbsp;&nbsp;       Read initial coords and vels from restart file
-  cut = 10.0;&nbsp;&nbsp;&nbsp;    Cutoff for non bonded interactions in Angstroms
-  ntr = 0;&nbsp;&nbsp;&nbsp;       No restrained atoms
-  ntc = 2;&nbsp;&nbsp;&nbsp;       SHAKE for constraining length of bonds involving Hydrogen atoms
-  ntf = 2;&nbsp;&nbsp;&nbsp;       Bond interactions involving H omitted
-  ntt = 3;&nbsp;&nbsp;&nbsp;       Constant temperature using Langevin dynamics
-  ig = -1;&nbsp;&nbsp;&nbsp;       Seed for pseudo-random number generator
-  ioutfm = 1;&nbsp;&nbsp;&nbsp;    Write trajectory in netcdf format
-  iwrap = 1;&nbsp;&nbsp;&nbsp;     Wrap coords into primary box
-  nstlim = 5000;&nbsp;&nbsp;&nbsp; Number of MD steps 
-  dt = 0.002;&nbsp;&nbsp;&nbsp;    Time step (in ps)
-  tempi = 0.0;&nbsp;&nbsp;&nbsp;   Initial temperature (0 K)
-  temp0 = 300.0;&nbsp;&nbsp;&nbsp; Final temperature (300 K)
-  irest = 0;&nbsp;&nbsp;&nbsp;     No restart from previous simulation
-  ntb = 1;&nbsp;&nbsp;&nbsp;       Periodic boundary conditions at constant volume
-  gamma_ln = 1.0;&nbsp;&nbsp;&nbsp;   Collision frequency for Langevin dynamics (in 1/ps)

In this particular example, the **heating** of the system is done in **2500 steps** (5ps) and is going **from 0K to 300K** (note that the number of steps has been reduced in this tutorial, for the sake of time). 


```python
# Import module
from biobb_amber.sander.sander_mdrun import sander_mdrun

# Create prop dict and inputs/outputs
output_heat_traj_path = 'sander.heat.netcdf'
output_heat_rst_path = 'sander.heat.rst'
output_heat_log_path = 'sander.heat.log'

prop = {
    "simulation_type" : "heat",
    "mdin" : { 
        'nstlim' : 2500, # Reducing the number of steps for the sake of time (5ps)
        'ntr' : 1,       # Overwritting restrain parameter
        'restraintmask' : '\"!:WAT,Cl-,Na+\"',      # Restraining solute
        'restraint_wt' : 10.0                       # With a force constant of 10 Kcal/mol*A2
    }
}

# Create and launch bb
sander_mdrun(input_top_path=output_ions_top_path,
            input_crd_path=output_min_rst_path,
            input_ref_path=output_min_rst_path,
            output_traj_path=output_heat_traj_path,
            output_rst_path=output_heat_rst_path,
            output_log_path=output_heat_log_path,
            properties=prop)
```

<a id="heatStep2"></a>
### Step 2: Checking results from the system warming up
Checking **system warming up** output. Plotting **temperature** along time during the **heating process**.


```python
# Import module
from biobb_amber.process.process_mdout import process_mdout

# Create prop dict and inputs/outputs
output_dat_heat_path = 'sander.md.temp.dat'

prop = {
    "terms" : ['TEMP']
}

# Create and launch bb
process_mdout(input_log_path=output_heat_log_path,
            output_dat_path=output_dat_heat_path,
            properties=prop)
```


```python
import plotly.graph_objs as go

with open(output_dat_heat_path, 'r') as energy_file:
    x, y = zip(*[
        (float(line.split()[0]), float(line.split()[1]))
        for line in energy_file
        if not line.startswith(("#", "@"))
        if float(line.split()[1]) < 1000
    ])

# Create a scatter plot
fig = go.Figure(data=go.Scatter(x=x, y=y, mode='lines'))

# Update layout
fig.update_layout(title="Heating process",
                  xaxis_title="Heating Step (ps)",
                  yaxis_title="Temperature (K)",
                  height=600)

# Show the figure (renderer changes for colab and jupyter)
rend = 'colab' if 'google.colab' in sys.modules else ''
fig.show(renderer=rend)
```

<img src='_static/md_setup/plot04.png'></img>

<a id="nvt"></a>
***
## Equilibrate the system (NVT)
Equilibrate the **protein system** in **NVT ensemble** (constant Number of particles, Volume and Temperature). Protein **heavy atoms** will be restrained using position restraining forces: movement is permitted, but only after overcoming a substantial energy penalty. The utility of position restraints is that they allow us to equilibrate our solvent around our protein, without the added variable of structural changes in the protein.

- [Step 1](#eqNVTStep1): Equilibrate the **protein system** with **NVT** ensemble.
- [Step 2](#eqNVTStep2): Checking **NVT Equilibration** results. Plotting **system temperature** by time during the **NVT equilibration** process.  
***
**Building Blocks** used:
 - [sander_mdrun](https://biobb-amber.readthedocs.io/en/latest/sander.html#module-sander.sander_mdrun) from **biobb_amber.sander.sander_mdrun**
 - [process_mdout](https://biobb-amber.readthedocs.io/en/latest/process.html#module-process.process_mdout) from **biobb_amber.process.process_mdout** 
***

<a id="eqNVTStep1"></a>
### Step 1: Equilibrating the system (NVT)
The **nvt** type of the **simulation_type property** contains the main default parameters to run a **system equilibration in NVT ensemble**:

-  imin = 0;&nbsp;&nbsp;&nbsp;      Run MD (no minimization)
-  ntx = 5;&nbsp;&nbsp;&nbsp;       Read initial coords and vels from restart file
-  cut = 10.0;&nbsp;&nbsp;&nbsp;    Cutoff for non bonded interactions in Angstroms
-  ntr = 0;&nbsp;&nbsp;&nbsp;       No restrained atoms
-  ntc = 2;&nbsp;&nbsp;&nbsp;       SHAKE for constraining length of bonds involving Hydrogen atoms
-  ntf = 2;&nbsp;&nbsp;&nbsp;       Bond interactions involving H omitted
-  ntt = 3;&nbsp;&nbsp;&nbsp;       Constant temperature using Langevin dynamics
-  ig = -1;&nbsp;&nbsp;&nbsp;       Seed for pseudo-random number generator
-  ioutfm = 1;&nbsp;&nbsp;&nbsp;    Write trajectory in netcdf format
-  iwrap = 1;&nbsp;&nbsp;&nbsp;     Wrap coords into primary box
-  nstlim = 5000;&nbsp;&nbsp;&nbsp; Number of MD steps 
-  dt = 0.002;&nbsp;&nbsp;&nbsp;    Time step (in ps)
-  irest = 1;&nbsp;&nbsp;&nbsp;     Restart previous simulation
-  ntb = 1;&nbsp;&nbsp;&nbsp;       Periodic boundary conditions at constant volume
-  gamma_ln = 5.0;&nbsp;&nbsp;&nbsp;   Collision frequency for Langevin dynamics (in 1/ps)

In this particular example, the **NVT equilibration** of the system is done in **500 steps** (note that the number of steps has been reduced in this tutorial, for the sake of time). 


```python
# Import module
from biobb_amber.sander.sander_mdrun import sander_mdrun

# Create prop dict and inputs/outputs
output_nvt_traj_path = 'sander.nvt.netcdf'
output_nvt_rst_path = 'sander.nvt.rst'
output_nvt_log_path = 'sander.nvt.log'

prop = {
    "simulation_type" : 'nvt',
    "mdin" : { 
        'nstlim' : 500, # Reducing the number of steps for the sake of time (1ps)
        'ntr' : 1,      # Overwritting restrain parameter
        'restraintmask' : '\"!:WAT,Cl-,Na+ & !@H=\"',      # Restraining solute heavy atoms
        'restraint_wt' : 5.0                               # With a force constant of 5 Kcal/mol*A2
    }
}

# Create and launch bb
sander_mdrun(input_top_path=output_ions_top_path,
            input_crd_path=output_heat_rst_path,
            input_ref_path=output_heat_rst_path,
            output_traj_path=output_nvt_traj_path,
            output_rst_path=output_nvt_rst_path,
            output_log_path=output_nvt_log_path,
            properties=prop)
```

<a id="eqNVTStep2"></a>
### Step 2: Checking NVT Equilibration results
Checking **NVT Equilibration** results. Plotting **system temperature** by time during the NVT equilibration process. 


```python
# Import module
from biobb_amber.process.process_mdout import process_mdout

# Create prop dict and inputs/outputs
output_dat_nvt_path = 'sander.md.nvt.temp.dat'

prop = {
    "terms" : ['TEMP']
}

# Create and launch bb
process_mdout(input_log_path=output_nvt_log_path,
            output_dat_path=output_dat_nvt_path,
            properties=prop)
```


```python
import plotly.graph_objs as go

with open(output_dat_nvt_path, 'r') as energy_file:
    x, y = zip(*[
        (float(line.split()[0]), float(line.split()[1]))
        for line in energy_file
        if not line.startswith(("#", "@"))
        if float(line.split()[1]) < 1000
    ])

# Create a scatter plot
fig = go.Figure(data=go.Scatter(x=x, y=y, mode='lines'))

# Update layout
fig.update_layout(title="NVT equilibration",
                  xaxis_title="Equilibration Step (ps)",
                  yaxis_title="Temperature (K)",
                  height=600)

# Show the figure (renderer changes for colab and jupyter)
rend = 'colab' if 'google.colab' in sys.modules else ''
fig.show(renderer=rend)
```

<img src='_static/md_setup/plot05.png'></img>

<a id="npt"></a>
***
## Equilibrate the system (NPT)
Equilibrate the **protein system** in **NPT ensemble** (constant Number of particles, Pressure and Temperature). Protein **heavy atoms** will be restrained using position restraining forces: movement is permitted, but only after overcoming a substantial energy penalty. The utility of position restraints is that they allow us to equilibrate our solvent around our protein, without the added variable of structural changes in the protein.

- [Step 1](#eqNPTStep1): Equilibrate the **protein system** with **NPT** ensemble.
- [Step 2](#eqNPTStep2): Checking **NPT Equilibration** results. Plotting **system pressure and density** by time during the **NVT equilibration** process.  
***
**Building Blocks** used:
 - [sander_mdrun](https://biobb-amber.readthedocs.io/en/latest/sander.html#module-sander.sander_mdrun) from **biobb_amber.sander.sander_mdrun**
 - [process_mdout](https://biobb-amber.readthedocs.io/en/latest/process.html#module-process.process_mdout) from **biobb_amber.process.process_mdout** 
***

<a id="eqNPTStep1"></a>
### Step 1: Equilibrating the system (NPT)
The **npt** type of the **simulation_type property** contains the main default parameters to run a **system equilibration in NPT ensemble**:

-  imin = 0;&nbsp;&nbsp;&nbsp;      Run MD (no minimization)
-  ntx = 5;&nbsp;&nbsp;&nbsp;       Read initial coords and vels from restart file
-  cut = 10.0;&nbsp;&nbsp;&nbsp;    Cutoff for non bonded interactions in Angstroms
-  ntr = 0;&nbsp;&nbsp;&nbsp;       No restrained atoms
-  ntc = 2;&nbsp;&nbsp;&nbsp;       SHAKE for constraining length of bonds involving Hydrogen atoms
-  ntf = 2;&nbsp;&nbsp;&nbsp;       Bond interactions involving H omitted
-  ntt = 3;&nbsp;&nbsp;&nbsp;       Constant temperature using Langevin dynamics
-  ig = -1;&nbsp;&nbsp;&nbsp;       Seed for pseudo-random number generator
-  ioutfm = 1;&nbsp;&nbsp;&nbsp;    Write trajectory in netcdf format
-  iwrap = 1;&nbsp;&nbsp;&nbsp;     Wrap coords into primary box
-  nstlim = 5000;&nbsp;&nbsp;&nbsp; Number of MD steps 
-  dt = 0.002;&nbsp;&nbsp;&nbsp;    Time step (in ps)
-  irest = 1;&nbsp;&nbsp;&nbsp;     Restart previous simulation
-  gamma_ln = 5.0;&nbsp;&nbsp;&nbsp;   Collision frequency for Langevin dynamics (in 1/ps)
-  pres0 = 1.0;&nbsp;&nbsp;&nbsp;   Reference pressure
-  ntp = 1;&nbsp;&nbsp;&nbsp;       Constant pressure dynamics: md with isotropic position scaling
-  taup = 2.0;&nbsp;&nbsp;&nbsp;    Pressure relaxation time (in ps)

In this particular example, the **NPT equilibration** of the system is done in **500 steps** (note that the number of steps has been reduced in this tutorial, for the sake of time). 


```python
# Import module
from biobb_amber.sander.sander_mdrun import sander_mdrun

# Create prop dict and inputs/outputs
output_npt_traj_path = 'sander.npt.netcdf'
output_npt_rst_path = 'sander.npt.rst'
output_npt_log_path = 'sander.npt.log'

prop = {
    "simulation_type" : 'npt',
    "mdin" : { 
        'nstlim' : 500, # Reducing the number of steps for the sake of time (1ps)
        'ntr' : 1,      # Overwritting restrain parameter
        'restraintmask' : '\"!:WAT,Cl-,Na+ & !@H=\"',      # Restraining solute heavy atoms
        'restraint_wt' : 2.5                               # With a force constant of 2.5 Kcal/mol*A2
    }
}

# Create and launch bb
sander_mdrun(input_top_path=output_ions_top_path,
            input_crd_path=output_nvt_rst_path,
            input_ref_path=output_nvt_rst_path,
            output_traj_path=output_npt_traj_path,
            output_rst_path=output_npt_rst_path,
            output_log_path=output_npt_log_path,
            properties=prop)
```

<a id="eqNPTStep2"></a>
### Step 2: Checking NPT Equilibration results
Checking **NPT Equilibration** results. Plotting **system pressure and density** by time during the **NPT equilibration** process. 


```python
# Import module
from biobb_amber.process.process_mdout import process_mdout

# Create prop dict and inputs/outputs
output_dat_npt_path = 'sander.md.npt.dat'

prop = {
    "terms" : ['PRES','DENSITY']
}

# Create and launch bb
process_mdout(input_log_path=output_npt_log_path,
            output_dat_path=output_dat_npt_path,
            properties=prop)
```


```python
# Read pressure and density data from file 
import plotly.graph_objs as go

# Read pressure and density data from file 
with open(output_dat_npt_path, 'r') as pd_file:
    x, y, z = zip(*[
        (float(line.split()[0]), float(line.split()[1]), float(line.split()[2]))
        for line in pd_file
        if not line.startswith(("#", "@"))
        if float(line.split()[1]) < 1000
    ])

# Create a scatter plot
trace1 = go.Scatter(
    x=x,y=y
)
trace2 = go.Scatter(
    x=x,y=z
)

fig = subplots.make_subplots(rows=1, cols=2, print_grid=False)

fig.append_trace(trace1, 1, 1)
fig.append_trace(trace2, 1, 2)

fig['layout']['xaxis1'].update(title='Time (ps)')
fig['layout']['xaxis2'].update(title='Time (ps)')
fig['layout']['yaxis1'].update(title='Pressure (bar)')
fig['layout']['yaxis2'].update(title='Density (Kg*m^-3)')

fig['layout'].update(title='Pressure and Density during NPT Equilibration')
fig['layout'].update(showlegend=False)
fig['layout'].update(height=500)

# Show the figure (renderer changes for colab and jupyter)
rend = 'colab' if 'google.colab' in sys.modules else ''
fig.show(renderer=rend)
```

<img src='_static/md_setup/plot06.png'></img>

<a id="free"></a>
***
## Free Molecular Dynamics Simulation
Upon completion of the **two equilibration phases (NVT and NPT)**, the system is now well-equilibrated at the desired temperature and pressure. The **position restraints** can now be released. The last step of the **protein** MD setup is a short, **free MD simulation**, to ensure the robustness of the system. 
- [Step 1](#mdStep1): Run short MD simulation of the **protein system**.
- [Step 2](#mdStep2): Checking results for the final step of the setup process, the **free MD run**. Plotting **Root Mean Square deviation (RMSd)** and **Radius of Gyration (Rgyr)** by time during the **free MD run** step.
***
**Building Blocks** used:
 - [sander_mdrun](https://biobb-amber.readthedocs.io/en/latest/sander.html#module-sander.sander_mdrun) from **biobb_amber.sander.sander_mdrun**
 - [process_mdout](https://biobb-amber.readthedocs.io/en/latest/process.html#module-process.process_mdout) from **biobb_amber.process.process_mdout** 
 - [cpptraj_rms](https://biobb-analysis.readthedocs.io/en/latest/ambertools.html#module-ambertools.cpptraj_rms) from **biobb_analysis.cpptraj.cpptraj_rms**
 - [cpptraj_rgyr](https://biobb-analysis.readthedocs.io/en/latest/ambertools.html#module-ambertools.cpptraj_rgyr) from **biobb_analysis.cpptraj.cpptraj_rgyr**
***

<a id="mdStep1"></a>
### Step 1: Creating portable binary run file to run a free MD simulation

The **free** type of the **simulation_type property** contains the main default parameters to run an **unrestrained MD simulation**:

-  imin = 0;&nbsp;&nbsp;&nbsp;      Run MD (no minimization)
-  ntx = 5;&nbsp;&nbsp;&nbsp;       Read initial coords and vels from restart file
-  cut = 10.0;&nbsp;&nbsp;&nbsp;    Cutoff for non bonded interactions in Angstroms
-  ntr = 0;&nbsp;&nbsp;&nbsp;       No restrained atoms
-  ntc = 2;&nbsp;&nbsp;&nbsp;       SHAKE for constraining length of bonds involving Hydrogen atoms
-  ntf = 2;&nbsp;&nbsp;&nbsp;       Bond interactions involving H omitted
-  ntt = 3;&nbsp;&nbsp;&nbsp;       Constant temperature using Langevin dynamics
-  ig = -1;&nbsp;&nbsp;&nbsp;       Seed for pseudo-random number generator
-  ioutfm = 1;&nbsp;&nbsp;&nbsp;    Write trajectory in netcdf format
-  iwrap = 1;&nbsp;&nbsp;&nbsp;     Wrap coords into primary box
-  nstlim = 5000;&nbsp;&nbsp;&nbsp; Number of MD steps 
-  dt = 0.002;&nbsp;&nbsp;&nbsp;    Time step (in ps)

In this particular example, a short, **5ps-length** simulation (2500 steps) is run, for the sake of time. 


```python
# Import module
from biobb_amber.sander.sander_mdrun import sander_mdrun

# Create prop dict and inputs/outputs
output_free_traj_path = 'sander.free.netcdf'
output_free_rst_path = 'sander.free.rst'
output_free_log_path = 'sander.free.log'

prop = {
    "simulation_type" : 'free',
    "mdin" : { 
        'nstlim' : 2500, # Reducing the number of steps for the sake of time (5ps)
        'ntwx' : 500  # Print coords to trajectory every 500 steps (1 ps)
    }
}

# Create and launch bb
sander_mdrun(input_top_path=output_ions_top_path,
            input_crd_path=output_npt_rst_path,
            output_traj_path=output_free_traj_path,
            output_rst_path=output_free_rst_path,
            output_log_path=output_free_log_path,
            properties=prop)
```

<a id="mdStep2"></a>
### Step 2: Checking free MD simulation results

Checking results for the final step of the setup process, the **free MD run**. Plotting **Root Mean Square deviation (RMSd)** and **Radius of Gyration (Rgyr)** by time during the **free MD run** step. **RMSd** against the **experimental structure** (input structure of the pipeline) and against the **minimized and equilibrated structure** (output structure of the NPT equilibration step).


```python
# cpptraj_rms: Computing Root Mean Square deviation to analyse structural stability 
#              RMSd against minimized and equilibrated snapshot (backbone atoms)   

# Import module
from biobb_analysis.ambertools.cpptraj_rms import cpptraj_rms

# Create prop dict and inputs/outputs
output_rms_first = pdbCode+'_rms_first.dat'

prop = {
    'mask': 'backbone',
    'reference': 'first'
}

# Create and launch bb
cpptraj_rms(input_top_path=output_ions_top_path,
            input_traj_path=output_free_traj_path,
            output_cpptraj_path=output_rms_first,
            properties=prop)
```


```python
# cpptraj_rms: Computing Root Mean Square deviation to analyse structural stability 
#              RMSd against experimental structure (backbone atoms)   

# Import module
from biobb_analysis.ambertools.cpptraj_rms import cpptraj_rms

# Create prop dict and inputs/outputs
output_rms_exp = pdbCode+'_rms_exp.dat'

prop = {
    'mask': 'backbone',
    'reference': 'experimental'
}

# Create and launch bb
cpptraj_rms(input_top_path=output_ions_top_path,
            input_traj_path=output_free_traj_path,
            output_cpptraj_path=output_rms_exp,
            input_exp_path=downloaded_pdb, 
            properties=prop)
```


```python
import plotly.graph_objs as go

# Read RMS vs first snapshot data from file 
with open(output_rms_first, 'r') as rms_first_file:
    x, y = zip(*[
        (float(line.split()[0]), float(line.split()[1]))
        for line in rms_first_file
        if not line.startswith(("#", "@"))
        if float(line.split()[1]) < 1000
    ])

# Read RMS vs experimental structure data from file 
with open(output_rms_exp, 'r') as rms_exp_file:
    x2, y2 = zip(*[
        (float(line.split()[0]), float(line.split()[1]))
        for line in rms_exp_file
        if not line.startswith(("#", "@"))
        if float(line.split()[1]) < 1000
    ])

# Create a scatter plot
fig = go.Figure(data=[go.Scatter(x=x, y=y, name = 'RMSd vs first'), go.Scatter(x=x, y=y2, name = 'RMSd vs exp')])

# Update layout
fig.update_layout(title="RMSd during free MD Simulation",
                  xaxis_title="Time (ps)",
                  yaxis_title="RMSd (Angstrom)",
                  height=600)

# Show the figure (renderer changes for colab and jupyter)
rend = 'colab' if 'google.colab' in sys.modules else ''
fig.show(renderer=rend)

```

<img src='_static/md_setup/plot07.png'></img>


```python
# cpptraj_rgyr: Computing Radius of Gyration to measure the protein compactness during the free MD simulation 

# Import module
from biobb_analysis.ambertools.cpptraj_rgyr import cpptraj_rgyr

# Create prop dict and inputs/outputs
output_rgyr = pdbCode+'_rgyr.dat'

prop = {
    'mask': 'backbone'
}

# Create and launch bb
cpptraj_rgyr(input_top_path=output_ions_top_path,
            input_traj_path=output_free_traj_path,
            output_cpptraj_path=output_rgyr,
            properties=prop)

```


```python
# Read Rgyr data from file 
import plotly.graph_objs as go

with open(output_rgyr, 'r') as rgyr_file:
    x, y = zip(*[
        (float(line.split()[0]), float(line.split()[1]))
        for line in rgyr_file
        if not line.startswith(("#", "@"))
        if float(line.split()[1]) < 1000
    ])

# Create a scatter plot
fig = go.Figure(data=go.Scatter(x=x, y=y, mode='lines'))

# Update layout
fig.update_layout(title="Radius of Gyration",
                  xaxis_title="Time (ps)",
                  yaxis_title="Rgyr (Angstrom)",
                  height=600)

# Show the figure (renderer changes for colab and jupyter)
rend = 'colab' if 'google.colab' in sys.modules else ''
fig.show(renderer=rend)
```

<img src='_static/md_setup/plot08.png'></img>

<a id="post"></a>
***
## Post-processing and Visualizing resulting 3D trajectory
Post-processing and Visualizing the **protein system** MD setup **resulting trajectory** using **NGL**
- [Step 1](#ppStep1): *Imaging* the resulting trajectory, **stripping out water molecules and ions** and **correcting periodicity issues**.
- [Step 2](#ppStep3): Visualizing the *imaged* trajectory using the *dry* structure as a **topology**. 
***
**Building Blocks** used:
 - [cpptraj_image](https://biobb-analysis.readthedocs.io/en/latest/ambertools.html#module-ambertools.cpptraj_image) from **biobb_analysis.cpptraj.cpptraj_image** 
***

<a id="ppStep1"></a>
### Step 1: *Imaging* the resulting trajectory.
Stripping out **water molecules and ions** and **correcting periodicity issues**  


```python
# cpptraj_image: "Imaging" the resulting trajectory
#                Removing water molecules and ions from the resulting structure

# Import module
from biobb_analysis.ambertools.cpptraj_image import cpptraj_image

# Create prop dict and inputs/outputs
output_imaged_traj = pdbCode+'_imaged_traj.trr'

prop = {
    'mask': 'solute',
    'format': 'trr'
}

# Create and launch bb
cpptraj_image(input_top_path=output_ions_top_path,
            input_traj_path=output_free_traj_path,
            output_cpptraj_path=output_imaged_traj,
            properties=prop)
```

<a id="ppStep2"></a>
### Step 2: Visualizing the generated dehydrated trajectory.
Using the **imaged trajectory** (output of the [Post-processing step 1](#ppStep1)) with the **dry structure** (output of the [Structure energy minimization - step 5](#minv)) as a topology.


```python
# Show trajectory
view = nglview.show_simpletraj(nglview.SimpletrajTrajectory(output_imaged_traj, output_ambpdb_path), gui=True)
view
```

<img src='_static/md_setup/traj.gif'></img>

<a id="output"></a>
## Output files

Important **Output files** generated:
 - structure.ions.pdb: **System structure** of the MD setup protocol. Structure generated during the MD setup and used in the MD simulation. With hydrogen atoms, solvent box and counterions.  
 - sander.free.netcdf: **Final trajectory** of the MD setup protocol.
 - sander.free.rst: **Final checkpoint file**, with information about the state of the simulation. It can be used to **restart** or **continue** a MD simulation.
 - structure.ions.parmtop: **Final topology** of the MD system in AMBER Parm7 format.

**Analysis** (MD setup check) output files generated:
 - 1aki_rms_first.dat: **Root Mean Square deviation (RMSd)** against **minimized and equilibrated structure** of the final **free MD run step**.
 - 1aki_rms_exp.dat: **Root Mean Square deviation (RMSd)** against **experimental structure** of the final **free MD run step**.
 - 1aki_rgyr.dat: **Radius of Gyration** of the final **free MD run step** of the **setup pipeline**.
 

***
<a id="questions"></a>

## Questions & Comments

Questions, issues, suggestions and comments are really welcome!

* GitHub issues:
    * [https://github.com/bioexcel/biobb](https://github.com/bioexcel/biobb)

* BioExcel forum:
    * [https://ask.bioexcel.eu/c/BioExcel-Building-Blocks-library](https://ask.bioexcel.eu/c/BioExcel-Building-Blocks-library)

