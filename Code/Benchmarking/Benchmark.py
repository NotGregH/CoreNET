
import subprocess

def run_benchmarks():

	filesList = ["1000genes.txt","2000genes.txt","3000genes.txt",
	"4000genes.txt","5000genes.txt","6000genes.txt",
	"7000genes.txt","8000genes.txt","9000genes.txt",
	"10000genes.txt","11000genes.txt","12000genes.txt",
	"13000genes.txt","14000genes.txt","15000genes.txt","All_Genes.txt"]

	for i in filesList:

		subprocess.call(["C:\WINDOWS\system32\WindowsPowerShell\\v1.0\powershell.exe", "python -m memory_profiler Network_Query.py -d Benchmark.db -o Test.SIF " + 
			"-g .\Test_Lists\Run2\\" + i + " | Out-File -Append  Run2_0Hop_noFilter.txt"], shell=True)

		subprocess.call(["C:\WINDOWS\system32\WindowsPowerShell\\v1.0\powershell.exe","python -m memory_profiler Network_Query.py -d Benchmark.db -o Test.SIF " + 
			"-j 1 -g .\Test_Lists\Run2\\" + i + " | Out-File -Append  Run2_1Hop_noFilter.txt"], shell = True)


		subprocess.call(["C:\WINDOWS\system32\WindowsPowerShell\\v1.0\powershell.exe","python -m memory_profiler Network_Query.py -d Benchmark.db -o Test.SIF " + 
			"-e .\Edge_lists\ppEdges.txt -g .\Test_Lists\Run2\\" + i + " | Out-File -Append  Run2_0Hop_Filter.txt"], shell=True)

		subprocess.call(["C:\WINDOWS\system32\WindowsPowerShell\\v1.0\powershell.exe","python -m memory_profiler Network_Query.py -d Benchmark.db -o Test.SIF " + 
			"-e .\Edge_lists\ppEdges.txt -j 1 -g .\Test_Lists\Run2\\" + i + " | Out-File -Append  Run2_1Hop_Filter.txt"], shell = True)


		subprocess.call(["C:\WINDOWS\system32\WindowsPowerShell\\v1.0\powershell.exe","python -m memory_profiler Network_Query.py -d Benchmark.db -o Test.SIF " + 
			"-g .\Test_Lists\Run3\\" + i + " | Out-File -Append  Run3_0Hop_noFilter.txt"], shell=True)

		subprocess.call(["C:\WINDOWS\system32\WindowsPowerShell\\v1.0\powershell.exe","python -m memory_profiler Network_Query.py -d Benchmark.db -o Test.SIF " + 
			"-j 1 -g .\Test_Lists\Run3\\" + i + " | Out-File -Append  Run3_1Hop_noFilter.txt"], shell=True)


		subprocess.call(["C:\WINDOWS\system32\WindowsPowerShell\\v1.0\powershell.exe","python -m memory_profiler Network_Query.py -d Benchmark.db -o Test.SIF " + 
			"-e .\Edge_lists\ppEdges.txt -g .\Test_Lists\Run3\\" + i + " | Out-File -Append  Run3_0Hop_Filter.txt"], shell=True)

		subprocess.call(["C:\WINDOWS\system32\WindowsPowerShell\\v1.0\powershell.exe","python -m memory_profiler Network_Query.py -d Benchmark.db -o Test.SIF " + 
			"-e .\Edge_lists\ppEdges.txt -j 1 -g .\Test_Lists\Run3\\" + i + " | Out-File -Append  Run3_1Hop_Filter.txt"], shell=True)


	return 

if __name__ == "__main__":
	run_benchmarks()