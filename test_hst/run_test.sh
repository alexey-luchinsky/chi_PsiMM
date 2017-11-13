nEv=1e5
echo "old"
../build/test.exe $1 $nEv ../src/my_old.dec $2_old
echo "new"
../build/test.exe $1 $nEv ../src/my_combines.dec $2_new

