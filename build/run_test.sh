rm old.txt new.txt
echo "old"
./test.exe $1 1e4 ../src/my_old.dec _old >> old.txt
echo "new"
./test.exe $1 1e4 ../src/my_combines.dec _new >>  new.txt
diff old.txt new.txt

