## check if files exist. It files are missing, produce error message. Otheriwse, run the first part
if [[ -f ./H5N1_VN1203_DE_Probes.txt && -f ./H5N1_VN1203_UNIVERSE_Probes.txt && -f KEGG_Pathway_Genes.txt ]]; then
	python3 Pham_project_1.py
else
	echo "Input files are missing!"
fi

## check if part 1 has completed by checking output file. If output file exists, run part 2
if [ -f ./pathway.txt ]; then
	echo "Part 1 completed!"
	python3 Pham_project_2.py
else
	echo "Part 1 did not complete!"
fi
