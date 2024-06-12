
import subprocess

def run_benchmarks():

	filesList = ["1000genes.txt","2000genes.txt","3000genes.txt",
	"4000genes.txt","5000genes.txt","6000genes.txt",
	"7000genes.txt","8000genes.txt","9000genes.txt",
	"10000genes.txt","11000genes.txt","12000genes.txt",
	"13000genes.txt","14000genes.txt","15000genes.txt"]

	for i in filesList:

		subprocess.call(["C:\WINDOWS\system32\WindowsPowerShell\\v1.0\powershell.exe", "python -m memory_profiler Network_Correlation.py " +
			"-d Benchmark.db -p 0.05 -e .\..\JoanDataV2.txt -o Test.SIF " + 
			"-g .\Test_Lists\\" + i + " | Out-File -Append  Correlation_Mem_noFilter.txt"], shell=True) 
		
		subprocess.call(["C:\WINDOWS\system32\WindowsPowerShell\\v1.0\powershell.exe", "python -m memory_profiler Network_Correlation.py " +
			"-d Benchmark.db -p 0.05 -e .\..\JoanDataV2.txt -o Test.SIF -l .\Edge_lists\ppEdges.txt " + 
			"-g .\Test_Lists\\" + i + " | Out-File -Append  Correlation_Mem_Filter.txt"], shell=True) 
	
	return 

if __name__ == "__main__":
	run_benchmarks()