
#!/usr/bin/python3
import subprocess
import os

#os.chdir("/Users/ibrahimahmed/October23/testProj/Result")


def run_edger(list1, path_):
    #os.chdir(path_)
    command = "/Library/Frameworks/R.framework/Resources/bin/Rscript"
    #script = "/Users/ibrahimahmed/Desktop/project2/script1.r"
    script = "/Users/ibrahimahmed/projects/GUI/script1.r"

    '''
     The separators are not recognized in 
     R script so commas are added for splitting text 
    '''

    # args = ["GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", "GSM1545538_purep53.txt", "GSM1545539_JMS8-2.txt",
            # "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt", "GSM1545542_JMS8-5.txt",  "GSM1545544_JMS9-P7c.txt",
            # "GSM1545545_JMS9-P8c.txt"]
    args1 = list1
    cmd = [command, script] + args1
    subprocess.check_output(cmd, universal_newlines=True)


# rscript()

def run_deseq2():#(list1, path_):
    # os.chdir(path_)
    command = "/Library/Frameworks/R.framework/Resources/bin/Rscript"
    script = "/Users/ibrahimahmed/projects/GUI/getDESeq2Data.r"
    #script = "/Users/ibrahimahmed/October23/testProj/pathway.r"

    '''
     The separators are not recognized in 
     R script so commas are added for splitting text 
    '''

    # args = ["GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", "GSM1545538_purep53.txt", "GSM1545539_JMS8-2.txt",
            # "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt", "GSM1545542_JMS8-5.txt",  "GSM1545544_JMS9-P7c.txt",
            # "GSM1545545_JMS9-P8c.txt"]

    #args1 = list1
    cmd = [command, script] #+ args1
    subprocess.check_output(cmd, universal_newlines=True)


def run_pathway(list1, path_):
    # os.chdir(path_)
    command = "/Library/Frameworks/R.framework/Resources/bin/Rscript"
    script = "/Users/ibrahimahmed/October23/testProj/pathway.r"

    '''
     The separators are not recognized in 
     R script so commas are added for splitting text 
    '''

    #args = ["GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", "GSM1545538_purep53.txt", "GSM1545539_JMS8-2.txt",
            #"GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt", "GSM1545542_JMS8-5.txt",  "GSM1545544_JMS9-P7c.txt",
            #"GSM1545545_JMS9-P8c.txt"]
    #args1 = list1
    cmd = [command, script] #+ args1
    subprocess.check_output(cmd, universal_newlines=True)

def run_emogea():
    command = "/Library/Frameworks/R.framework/Resources/bin/Rscript"
    script = "/Users/ibrahimahmed/Desktop/emogea/GerlinskiDiffExpression1401.R"

    cmd = [command, script] #+ args1
    subprocess.check_output(cmd, universal_newlines=True)

if __name__=="__main__":
    run_deseq2()