../../../scripts/gen_particles.py -n 1000 particles.json
../../../scripts/gen_cloth.py -n 50 -ks 5000 -kd 50 -s 1 -f '{"o":[0,0,0],"x":[1,0,0],"y":[0,0,-1],"z":[0,1,0]}' cloth01.json
../../../scripts/gen_cloth.py -n 50 -ks 5000 -kd 50 -s 1 -p -f '{"o":[0,0,0],"x":[1,0,0],"y":[0,0,-1],"z":[0,1,0]}' cloth02.json
../../../scripts/gen_cloth.py -n 50 -ks 5000 -kd 50 -s 1 cloth03.json
