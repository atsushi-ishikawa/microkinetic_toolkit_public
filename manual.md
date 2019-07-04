### Methane

* input file = "input.txt"

#### 1. Run vasp
edit `run.sh` to change input file and then
```bash
qsub run.sh
```

#### 2. Make species file
```bash
python count_species.py input.txt
```

#### 3. Make pre-exponential file
```bash
python pre_exponential.py input.txt
```

#### 4. Make entropy file
```bash
python entropy.py input.txt
```

#### 5. Make MATLAB ODE file
```bash
python make_rate_equation.py input.txt
```

#### 6. Run MATLAB
```bash
cp pre_exp.txt met001ode.m deltaE.txt deltaS.txt barrier.txt species.txt  MATLAB_dir
```
* `coverage.txt` and `rateconst.txt` are generated

#### 7. Add rate to reactionfile
```bash
python rate_for_graph.py [input.txt] coverage.txt rateconst.txt variables.txt
```
* `input_val.txt` is generated

#### 8. Draw graph
```bash
python draw_chemical_graph.py input_val.txt coverage.txt
```

### Adsorbate specification
* cation and anion
CH3^+1 or CH3^-1

* multiple adsorbate on one surface
H_atop.x1y1,H_atop.x2y2 --> H2 + 2\*surf

* flip, side, high
CH4-FLIP_atop.x1y1 --> CH4 + surf
CH4-SIDE_atop.x1y1 --> CH4 + surf
CH4-HIGH_atop.x1y1 --> CH4 + surf
