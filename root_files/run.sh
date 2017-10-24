echo "./run nEv ff chi"
echo "The number of parameters is $#"
nEv=$1
ff=$2
chi=$3
../build/chiC_all.exe ${chi}0 $nEv ../src/my_$ff.dec _$ff
../build/chiC_all.exe ${chi}1 $nEv ../src/my_$ff.dec _$ff
../build/chiC_all.exe ${chi}2 $nEv ../src/my_$ff.dec _$ff
