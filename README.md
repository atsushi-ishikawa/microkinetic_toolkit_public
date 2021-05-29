### DFT計算を行い反応エネルギーを算出するとともにMATLABコードを生成する計算用スクリプト(Ref. ACS. Catal)

* インプットファイル: input.txt"
素反応を書く。詳しい書き方は下に。

* ステップ
#### 1. VASP計算を行う
`run.sh`を編集
```bash
qsub or pjsub run.sh
```

#### 2. 化学種をカウントする
```bash
python count_species.py input.txt
```

#### 3. pere-exponentialを作る
```bash
python pre_exponential.py input.txt
```

#### 4. entropyを計算
```bash
python entropy.py input.txt
```

#### 5. MATLAB ODEファイルと作る
```bash
python make_rate_equation.py input.txt
```

#### 6. MATLAB計算を行う
```bash
cp pre_exp.txt met001ode.m deltaE.txt deltaS.txt barrier.txt species.txt  MATLAB_dir
```
* `coverage.txt` and `rateconst.txt` are generated

#### 7. Add rate to reactionfile
```bash
python rate_for_graph.py [input.txt] coverage.txt rateconst.txt variables.txt
```
* `input_val.txt` is generated
or 
```bash
python rate_for_graph2.py [input.txt] nodes.txt rateconst.txt variables.txt
```
* `edges.txt` is generated

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
CH4-FLIP_atop.x1y1  --> CH4 + surf
CH4-SIDEx_atop.x1y1 --> CH4 + surf
CH4-SIDEy_atop.x1y1 --> CH4 + surf
CH4-TILT_atop.x1y1  --> CH4 + surf
CH4-HIGH_atop.x1y1  --> CH4 + surf
CH4-AIR_atop.x1y1   --> CH4 + surf
