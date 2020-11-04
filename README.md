# simraStudies

mkdir -p ~/kode
mkdir -p ~/results/simra/HuntHill
cd ~/kode
git clone https://github.com/Zetison/colormaps
git clone https://github.com/Zetison/paraUtils

# Run SIMRA
```
python3 generate_terrain.py

gridgen

gfortran ic_bc.f90 -o ic_bc

./ic_bc simra.in

simra simra.in
```
