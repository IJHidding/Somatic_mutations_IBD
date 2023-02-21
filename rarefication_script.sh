#Gene=$1
for file in $(cat biopsies_files.txt); 
do 
	grep -v '^#' $file > stratification/variable_$file
	python3 recalculate.py stratification/variable_$file
done
