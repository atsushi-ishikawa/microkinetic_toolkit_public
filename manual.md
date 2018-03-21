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

#### 4. Make MATLAB ODE file
```bash
python make_rate_equation.py input.txt
```
