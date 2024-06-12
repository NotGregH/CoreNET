
import subprocess

def run_benchmarks():


	subprocess.call(["C:\WINDOWS\system32\WindowsPowerShell\\v1.0\powershell.exe", "python .\CrossSpeciesTool.py " \
		"-e1 ..\ArabidopsisData\\filteredgenes_NormExpMas5.tab -d1 ..\ArabCross.sqlite -oM ..\orthoFiles\OrthoMCL.groups.txt -p 0.05 -s2 " \
		"..\MaizeData\MaizeMultiNet.sif -e2 ..\MaizeData\uniqueGenes_filteredNormalized_exp1.txt -d B -o MaizeExp1AraMCL.sif AraMaizeExp1MCL.sif " \
		"| Out-file -Append CrossToolTime.txt"], shell=True) 
		

	subprocess.call(["C:\WINDOWS\system32\WindowsPowerShell\\v1.0\powershell.exe", "python .\CrossSpeciesTool.py " \
		"-e1 ..\ArabidopsisData\\filteredgenes_NormExpMas5.tab -d1 ..\ArabCross.sqlite -oM ..\orthoFiles\OrthoMCL.groups.txt -p 0.05 -s2 " \
		"..\MaizeData\MaizeMultiNet.sif -e2 ..\MaizeData\Exp2\genes_filteredNormalized_exp2.txt -d B -o MaizeExp2AraMCL.sif AraMaizeExp2MCL.sif " \
		"| Out-file -Append CrossToolTime.txt"], shell=True) 
		

	subprocess.call(["C:\WINDOWS\system32\WindowsPowerShell\\v1.0\powershell.exe", "python .\CrossSpeciesTool.py " \
		"-e1 ..\ArabidopsisData\\filteredgenes_NormExpMas5.tab -d1 ..\ArabCross.sqlite -oM ..\orthoFiles\OrthoMCL.groups.txt -p 0.05 -s2 " \
		"..\MaizeData\MaizeMultiNet.sif -e2 ..\MaizeData\Exp4\genes_Significant.probes.E4.N.txt -d B -o MaizeExp4AraMCL.sif AraMaizeExp4MCL.sif " \
		"| Out-file -Append CrossToolTime.txt"], shell=True) 
		

	subprocess.call(["C:\WINDOWS\system32\WindowsPowerShell\\v1.0\powershell.exe", "python .\CrossSpeciesTool.py " \
		"-e1 ..\ArabidopsisData\\filteredgenes_NormExpMas5.tab -d1 ..\ArabCross.sqlite -oM ..\orthoFiles\OrthoMCL.groups.txt -p 0.05 -s2 " \
		"..\MaizeData\MaizeMultiNet.sif -e2 ..\MaizeData\Exp5\\nitrogen_only_2554UniqueGenesExprs.txt -d B -o MaizeExp5AraMCL.sif AraMaizeExp5MCL.sif " \
		"| Out-file -Append CrossToolTime.txt"], shell=True) 
	
	return 

if __name__ == "__main__":
	run_benchmarks()
