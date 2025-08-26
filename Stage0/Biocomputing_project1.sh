echo "===== Project 1 ====="
#1. Print your name
$echo "Chama"
#2. Create a folder titled your name
$mkdir chama
#3. Create another new directory titled biocomputing and change to that directory with one line of command
$mkdir biocomputing && cd biocomputing
#4. Download these 3 files:
$wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.fna
$wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk
$wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk
#5. Move the .fna file to the folder titled your name
$mv *.fna ~/chama
#6. Delete the duplicate gbk file
$cd .. && cd biocomputing
$rm wildtype.gbk.1
#7. Confirm if the .fna file is mutant or wild type (tatatata vs tata)
$cd .. && cd chama
$grep -o "tatatata" wildtype.fna | wc -l
$grep -o "tata" wildtype.fna | wc -l
#8. If mutant, print all matching lines into a new file
$grep -i "tata" wildtype.fna > mutant_match.fna
#9. Count number of lines (excluding header) in the .gbk file
$cd .. && cd biocomputing
$grep -A3291 "^ORIGIN" wildtype.gbk | tail -n +2 | wc -l
#10. Print the sequence length of the .gbk file. (Use the LOCUS tag in the first line)
$grep "^LOCUS" wildtype.gbk | awk '{print $3}'
#11. Print the source organism of the .gbk file. (Use the SOURCE tag in the first line)
$grep "^SOURCE" wildtype.gbk | awk '{print $2,  $3}'
#12. List all the gene names of the .gbk file. Hint {grep '/gene='}
$grep "/gene=" wildtype.gbk
#13. Clear your terminal space and print all commands used today
$clear && history
#14. List the files in the two folders and share a screenshot of your terminal
$cd
$ls biocomputing/ && ls chama/
