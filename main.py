# The main script to manage the subworkflows for RNAseq data analysis
import yaml
import os
import time
from tkinter import *

#root = Tk()  

#def myClick():
    #datapath = "Data folder" + e.get()
    #myLabel = Label(root, text= datapath)
    #myLabel.pack()

#myButton = Button(root,text="Lexogen QuantSeq 3' Fwd",padx=25,command=myClick,fg="blue",bg="green")
#myButton.pack()  
##myButton.grid(row=0,column=5)

#e = Entry(root, width=100)
#e.pack()
#e.insert(0,"Path to data folder: ")
#root.mainloop()

with open('configs/config_main.yaml') as yamlfile:
    config = yaml.full_load(yamlfile)

#subprocess.run('source activate environment-name && "enter command here" && source deactivate', shell=True)

# Parameters to control the workflow
project = config["PROJECT"]

### Verify the reference and the annotation for it
reference = config["GENOME"]
print("--->>Is the reference genome correct? \n", reference)

genome_annotation = config["ANNOTATION"]
print("--->>Is the reference genome annotation correct? \n", genome_annotation)

# Start the workflow
print("Start rnaseq workflow on project: " + project)

## write the running time in a log file
file_log_time = open("logs/log_running_time.txt", "a+")
file_log_time.write("\nProject name: " + project + "\n")
file_log_time.write("Start time: " + time.ctime() + "\n")

def spend_time(start_time, end_time):
    seconds = end_time - start_time
    hours = seconds // 3600
    seconds %= 3600
    minutes = seconds // 60
    seconds %= 60

    return "%d:%02d:%02d" % (hours, minutes, seconds)

print("Start Quality Control!")
start_time = time.time()
#os.system("/opt/software/miniconda3_py38/bin/snakemake  -p --latency-wait 36000  --cluster \"qsub -V -l nodes=1:ppn=4,mem=40gb\" --jobs 10 -s workflow/quality_control.rules 2>&1 | tee logs/log_quality_control.txt")
#os.system("snakemake  -p --latency-wait 36000  --cluster \"qsub -l nodes=1:ppn=4,mem=40gb\" --jobs 10 -s workflow/quality_control.rules --dag | dot -Tsvg > logs/quality_control.svg")
os.system("snakemake  -p --latency-wait 36000  --jobs 4 -s workflow/quality_control.rules 2>&1 | tee logs/log_quality_control.txt")
end_time = time.time()
file_log_time.write("Time of running QC: " + spend_time(start_time, end_time) + "\n")


print("Start Trimming!")
start_time = time.time()
#os.system("/opt/software/miniconda3_py38/bin/snakemake -p --latency-wait 36000  --cluster \"qsub -V -l nodes=1:ppn=4,mem=40gb\" --jobs 10 -s workflow/trim.rules 2>&1 | tee logs/log_trim.txt")
os.system("snakemake -p --latency-wait 360  --rerun-incomplete --jobs 4 -s workflow/trim.rules 2>&1 | tee logs/log_trim.txt")
#os.system("snakemake -p --latency-wait 36000  --cluster \"qsub -l nodes=1:ppn=4,mem=40gb\" --jobs 10 -s workflow/trim.rules --dag | dot -Tsvg > logs/trim.svg")

end_time = time.time()
file_log_time.write("Time of running trimming:" + spend_time(start_time, end_time) + "\n")
print("Trimming is done!")

# Do alignment, QC on the alignment and feature extraction
print("Start Aligning and Feature extraction!")
start_time = time.time()
#os.system("snakemake -p --latency-wait 36000 --cluster \"qsub -V -l nodes=1:ppn=4,mem=80gb\" --jobs 10 -s workflow/alignment.rules 2>&1 | tee logs/log_align_count_genome.txt")
os.system("snakemake -p --latency-wait 36000  --jobs 3 -s workflow/alignment.rules 2>&1 | tee logs/log_align_count_genome.txt")
#os.system("snakemake -p --latency-wait 36000 --cluster \"qsub -l nodes=1:ppn=4,mem=80gb\" --jobs 10 -s workflow/alignment.rules --dag |dot -Tsvg > logs/align.svg")
end_time = time.time()
#file_log_time.write("Time of running genome alignment:" + spend_time(start_time, end_time) + "\n")

# Do differential gene expression
print("Start doing DEA!")
start_time = time.time()
#os.system("/opt/software/miniconda3_py38/bin/snakemake -p --latency-wait 36000 --cluster \"qsub -V -l nodes=1:ppn=4,mem=40gb\"  --jobs 10 -s workflow/dea.rules 2>&1 | tee logs/log_dea.txt")
#os.system("snakemake -p --latency-wait 36000 --cluster \"qsub -l nodes=1:ppn=4,mem=40gb\"  --jobs 10 -s workflow/dea.rules --dag | dot -Tsvg > logs/dea.svg")
tart_time = time.time()
os.system("snakemake -p --latency-wait 36000  -s workflow/dea.rules 2>&1 | tee logs/log_dea.txt")
end_time = time.time()
file_log_time.write("Time of running DEA based:" + spend_time(start_time, end_time) + "\n")
print("DEA is done!")


file_log_time.write("Finish time: " + time.ctime() + "\n")
file_log_time.close()
