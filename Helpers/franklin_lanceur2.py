#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Author: Luc Marchand
#This script will be run every hour to check if new LowPass data are present in 2 folders of /staging/hiseq_raw and if so, start the processing of the run.
#Meant to be used with Franklin13.py and above.  Work folders are created here instead of inside franklin13.py (and above).

#import urllib.request, json
import os
import glob
import sys
import subprocess
import time

def sendemail(emailto,run,message,file):
    if file == "":
        command = "echo \"" + message + "\" | mutt -s \"LowPass run " + run + "\" -e 'my_hdr From:LowPass Bioinfo <from@address.com>' -e 'set sendmail=/usr/sbin/ssmtp' -- " + emailto
    else:
        command = "echo \"" + message + "\" | mutt -s \"LowPass run " + run + "\" -a " + file + " -e 'my_hdr From:LowPass Bioinfo <from@address.com>' -e 'set sendmail=/usr/sbin/ssmtp' -- " + emailto
    
    tries = 0
    rcode = 99
    while rcode != 0 and tries < 5:
        tries = tries + 1
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
        process.wait()
        rcode = process.returncode
        if rcode != 0:
            print("Error code:")
            print(rcode)
            time.sleep(10)

emailto = "alexandre.dionne-laporte.hsj@ssss.gouv.qc.ca luc.marchand.hsj@ssss.gouv.qc.ca quang-hien.le.hsj@ssss.gouv.qc.ca dan.spiegelman.hsj@ssss.gouv.qc.ca emilie.lacroix.hsj@ssss.gouv.qc.ca"
file = ""

#First location A00516
dirs = os.listdir("/staging/hiseq_raw/A00516/")
for x in dirs:
    if os.path.isfile("/staging/hiseq_raw/A00516/" + x + "/CopyComplete.txt"):
        if not (os.path.isfile("/staging/hiseq_raw/A00516/" + x + "/processed.txt") or os.path.isfile("/staging/hiseq_raw/A00516/" + x + "/analyzing.txt") or os.path.isfile("/staging/hiseq_raw/A00516/" + x + "/failed.txt")):
            if len(glob.glob("/staging/hiseq_raw/A00516/" + x + "/LowPass*.csv")) == 1:
                print("New low pass run to analyze found: " + x )
                runname = x
                os.system("mkdir /staging2/dragen/" + runname)
                os.system("mkdir /staging2/dragen/" + runname + "/Run_logs")
                message = "Départ d'analyse LowPass"
                sendemail(emailto,runname,message,file)
                #os.system("echo \"Départ d'analyse LowPass\" | mutt -s \"Run Lowpass " + runname + "\" -e 'my_hdr From:LowPass Bioinfo <from@address.com>' -e 'set sendmail=/usr/sbin/ssmtp' -- luc.marchand.hsj@ssss.gouv.qc.ca")
                os.system("touch /staging/hiseq_raw/A00516/" + x + "/analyzing.txt")
                #command = "cd $HOME/" + runname + " && script -c \"python /staging2/soft/CHUSJ-utils/franklin13.py " + runname + "\" -f progress.txt"
                command = "cd /staging2/dragen/" + runname + "/Run_logs && python /staging2/soft/CHUSJ-utils/franklin15.py " + runname + " > progress.txt"

                #command = "cd /staging2/dragen/" + runname + "/Run_logs && script -c \"python /home/marluc00/franklintestlanceur.py " + runname + "\" -f progress.txt"
                #command = "cd /staging2/dragen/" + runname + "/Run_logs && python /home/marluc00/franklintestlanceur.py " + runname + " > progress.txt"

                process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
                process.wait()
                print(process.returncode)
                if process.returncode != 0:
                    print("L'analyse de la run " + runname + " a eu un problème.  Vérifier les logs.")
                    os.system("touch /staging/hiseq_raw/A00516/" + x + "/failed.txt")
                    os.system("rm /staging/hiseq_raw/A00516/" + x + "/analyzing.txt")
                    message = "Erreur dans la run, vérifier les logs"
                    sendemail(emailto,runname,message,file)
                    #os.system("echo \"Erreur dans la run, vérifier les logs\" | mutt -s \"Run Lowpass " + runname + "\" -e 'my_hdr From:LowPass Bioinfo <from@address.com>' -e 'set sendmail=/usr/sbin/ssmtp' -- luc.marchand.hsj@ssss.gouv.qc.ca")
                else:
                    os.system("touch /staging/hiseq_raw/A00516/" + x + "/processed.txt")
                    os.system("rm /staging/hiseq_raw/A00516/" + x + "/analyzing.txt")
                    message = "Analyse LowPass terminée"
                    sendemail(emailto,runname,message,file)
                    #os.system("echo \"Analyse LowPass terminée\" | mutt -s \"Run Lowpass " + runname + "\" -e 'my_hdr From:LowPass Bioinfo <from@address.com>' -e 'set sendmail=/usr/sbin/ssmtp' -- luc.marchand.hsj@ssss.gouv.qc.ca")
                #os.system("cd $HOME/" + runname + " && script -c \"python /staging2/soft/CHUSJ-utils/franklin12.py " + runname + "\" -f progress.txt")
                print("Fin de la pipeline")

            elif len(glob.glob("/staging/hiseq_raw/A00516/" + x + "/LowPass*.csv")) > 1:
                message = "Il y a plus qu'une SampleSheet dans la run"
                sendemail(emailto,runname,message,file)
                #os.system("echo \"Il y a plus qu'une SampleSheet dans la run\" | mutt -s \"Run Lowpass " + x + "\" -e 'my_hdr From:LowPass Bioinfo <from@address.com>' -- luc.marchand.hsj@ssss.gouv.qc.ca")

#Second location LH00336
dirs = os.listdir("/staging/hiseq_raw/LH00336/")
for x in dirs:
    if os.path.isfile("/staging/hiseq_raw/LH00336/" + x + "/CopyComplete.txt"):
        if not (os.path.isfile("/staging/hiseq_raw/LH00336/" + x + "/processed.txt") or os.path.isfile("/staging/hiseq_raw/LH00336/" + x + "/analyzing.txt") or os.path.isfile("/staging/hiseq_raw/LH00336/" + x + "/failed.txt")):
            if os.path.isfile("/staging/hiseq_raw/LH00336/" + x + "/SampleSheet.csv"):
                file1 = open("/staging/hiseq_raw/LH00336/" + x + "/SampleSheet.csv")
                for y in file1:
                    y = y.rstrip()
                    if "RunName" in y and "LowPass" in y:
                    #if len(glob.glob("/staging/hiseq_raw/LH00336/" + x + "/LowPass*.csv")) == 1:
                        if os.path.isfile("/staging/hiseq_raw/LH00336/" + x + "/Analysis/1/CopyComplete.txt"):
                            print("New low pass run to analyze found: " + x )
                            runname = x
                            os.system("mkdir /staging2/dragen/" + runname)
                            os.system("mkdir /staging2/dragen/" + runname + "/Run_logs")
                            message = "Départ d'analyse LowPass"
                            sendemail(emailto,runname,message,file)
                            #os.system("echo \"Départ d'analyse LowPass\" | mutt -s \"Run Lowpass " + runname + "\" -e 'my_hdr From:LowPass Bioinfo <from@address.com>' -e 'set sendmail=/usr/sbin/ssmtp' -- luc.marchand.hsj@ssss.gouv.qc.ca")
                            os.system("touch /staging/hiseq_raw/LH00336/" + x + "/analyzing.txt")
                            command = "cd /staging2/dragen/" + runname + "/Run_logs && python /staging2/soft/CHUSJ-utils/franklin15.py " + runname + " > progress.txt"

                            process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
                            process.wait()
                            print(process.returncode)
                            if process.returncode != 0:
                                print("L'analyse de la run " + runname + " a eu un problème.  Vérifier les logs.")
                                os.system("touch /staging/hiseq_raw/LH00336/" + x + "/failed.txt")
                                os.system("rm /staging/hiseq_raw/LH00336/" + x + "/analyzing.txt")
                                message = "Erreur dans la run, vérifier les logs"
                                sendemail(emailto,runname,message,file)
                                #os.system("echo \"Erreur dans la run, vérifier les logs\" | mutt -s \"Run Lowpass " + runname + "\" -e 'my_hdr From:LowPass Bioinfo <from@address.com>' -e 'set sendmail=/usr/sbin/ssmtp' -- luc.marchand.hsj@ssss.gouv.qc.ca")
                            else:
                                os.system("touch /staging/hiseq_raw/LH00336/" + x + "/processed.txt")
                                os.system("rm /staging/hiseq_raw/LH00336/" + x + "/analyzing.txt")
                                message = "Analyse LowPass terminée"
                                sendemail(emailto,runname,message,file)
                                #os.system("echo \"Analyse LowPass terminée\" | mutt -s \"Run Lowpass " + runname + "\" -e 'my_hdr From:LowPass Bioinfo <from@address.com>' -e 'set sendmail=/usr/sbin/ssmtp' -- luc.marchand.hsj@ssss.gouv.qc.ca")
                            #os.system("cd $HOME/" + runname + " && script -c \"python /staging2/soft/CHUSJ-utils/franklin12.py " + runname + "\" -f progress.txt")
                            print("Fin de la pipeline")
            else:
                print("Aucune SampleSheet trouvée pour la run " + x)
                os.system("touch /staging/hiseq_raw/LH00336/" + x + "/failed.txt")

#Third location LH00207R (Genome Quebec NovaseqX)
dirs = os.listdir("/staging/hiseq_raw/LH00207R/")
for x in dirs:
    if os.path.isfile("/staging/hiseq_raw/LH00207R/" + x + "/CopyComplete.txt"):
        if not (os.path.isfile("/staging/hiseq_raw/LH00207R/" + x + "/processed.txt") or os.path.isfile("/staging/hiseq_raw/LH00207R/" + x + "/analyzing.txt") or os.path.isfile("/staging/hiseq_raw/LH00207R/" + x + "/failed.txt")):
            if os.path.isfile("/staging/hiseq_raw/LH00207R/" + x + "/SampleSheet.csv"):
                file1 = open("/staging/hiseq_raw/LH00207R/" + x + "/SampleSheet.csv")
                for y in file1:
                    y = y.rstrip()
                    if "RunName" in y and "LowPass" in y:
                    #if len(glob.glob("/staging/hiseq_raw/LH00207R/" + x + "/LowPass*.csv")) == 1:
                        if os.path.isfile("/staging/hiseq_raw/LH00207R/" + x + "/Analysis/1/CopyComplete.txt"):
                            print("New low pass run to analyze found: " + x )
                            runname = x
                            os.system("mkdir /staging2/dragen/" + runname)
                            os.system("mkdir /staging2/dragen/" + runname + "/Run_logs")
                            message = "Départ d'analyse LowPass"
                            sendemail(emailto,runname,message,file)
                            #os.system("echo \"Départ d'analyse LowPass\" | mutt -s \"Run Lowpass " + runname + "\" -e 'my_hdr From:LowPass Bioinfo <from@address.com>' -e 'set sendmail=/usr/sbin/ssmtp' -- luc.marchand.hsj@ssss.gouv.qc.ca")
                            os.system("touch /staging/hiseq_raw/LH00207R/" + x + "/analyzing.txt")
                            command = "cd /staging2/dragen/" + runname + "/Run_logs && python /staging2/soft/CHUSJ-utils/franklin15.py " + runname + " > progress.txt"

                            process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
                            process.wait()
                            print(process.returncode)
                            if process.returncode != 0:
                                print("L'analyse de la run " + runname + " a eu un problème.  Vérifier les logs.")
                                os.system("touch /staging/hiseq_raw/LH00207R/" + x + "/failed.txt")
                                os.system("rm /staging/hiseq_raw/LH00207R/" + x + "/analyzing.txt")
                                message = "Erreur dans la run, vérifier les logs"
                                sendemail(emailto,runname,message,file)
                                #os.system("echo \"Erreur dans la run, vérifier les logs\" | mutt -s \"Run Lowpass " + runname + "\" -e 'my_hdr From:LowPass Bioinfo <from@address.com>' -e 'set sendmail=/usr/sbin/ssmtp' -- luc.marchand.hsj@ssss.gouv.qc.ca")
                            else:
                                os.system("touch /staging/hiseq_raw/LH00207R/" + x + "/processed.txt")
                                os.system("rm /staging/hiseq_raw/LH00207R/" + x + "/analyzing.txt")
                                message = "Analyse LowPass terminée"
                                sendemail(emailto,runname,message,file)
                                #os.system("echo \"Analyse LowPass terminée\" | mutt -s \"Run Lowpass " + runname + "\" -e 'my_hdr From:LowPass Bioinfo <from@address.com>' -e 'set sendmail=/usr/sbin/ssmtp' -- luc.marchand.hsj@ssss.gouv.qc.ca")
                            #os.system("cd $HOME/" + runname + " && script -c \"python /staging2/soft/CHUSJ-utils/franklin12.py " + runname + "\" -f progress.txt")
                            print("Fin de la pipeline")
            else:
                print("Aucune SampleSheet trouvée pour la run " + x)
                os.system("touch /staging/hiseq_raw/LH00207R/" + x + "/failed.txt")


#Fourth location A00977
dirs = os.listdir("/staging/hiseq_raw/A00977/")
for x in dirs:
    if os.path.isfile("/staging/hiseq_raw/A00977/" + x + "/CopyComplete.txt"):
        if not (os.path.isfile("/staging/hiseq_raw/A00977/" + x + "/processed.txt") or os.path.isfile("/staging/hiseq_raw/A00977/" + x + "/analyzing.txt") or os.path.isfile("/staging/hiseq_raw/A00977/" + x + "/failed.txt")):
            if len(glob.glob("/staging/hiseq_raw/A00977/" + x + "/LowPass*.csv")) == 1:
                print("New low pass run to analyze found: " + x )
                runname = x
                os.system("mkdir /staging2/dragen/" + runname)
                os.system("mkdir /staging2/dragen/" + runname + "/Run_logs")
                message = "Départ d'analyse LowPass"
                sendemail(emailto,runname,message,file)
                #os.system("echo \"Départ d'analyse LowPass\" | mutt -s \"Run Lowpass " + runname + "\" -e 'my_hdr From:LowPass Bioinfo <from@address.com>' -e 'set sendmail=/usr/sbin/ssmtp' -- luc.marchand.hsj@ssss.gouv.qc.ca")
                os.system("touch /staging/hiseq_raw/A00977/" + x + "/analyzing.txt")
                #command = "cd $HOME/" + runname + " && script -c \"python /staging2/soft/CHUSJ-utils/franklin13.py " + runname + "\" -f progress.txt"
                command = "cd /staging2/dragen/" + runname + "/Run_logs && python /staging2/soft/CHUSJ-utils/franklin15.py " + runname + " > progress.txt"

                #command = "cd /staging2/dragen/" + runname + "/Run_logs && script -c \"python /home/marluc00/franklintestlanceur.py " + runname + "\" -f progress.txt"
                #command = "cd /staging2/dragen/" + runname + "/Run_logs && python /home/marluc00/franklintestlanceur.py " + runname + " > progress.txt"

                process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
                process.wait()
                print(process.returncode)
                if process.returncode != 0:
                    print("L'analyse de la run " + runname + " a eu un problème.  Vérifier les logs.")
                    os.system("touch /staging/hiseq_raw/A00977/" + x + "/failed.txt")
                    os.system("rm /staging/hiseq_raw/A00977/" + x + "/analyzing.txt")
                    message = "Erreur dans la run, vérifier les logs"
                    sendemail(emailto,runname,message,file)
                    #os.system("echo \"Erreur dans la run, vérifier les logs\" | mutt -s \"Run Lowpass " + runname + "\" -e 'my_hdr From:LowPass Bioinfo <from@address.com>' -e 'set sendmail=/usr/sbin/ssmtp' -- luc.marchand.hsj@ssss.gouv.qc.ca")
                else:
                    os.system("touch /staging/hiseq_raw/A00977/" + x + "/processed.txt")
                    os.system("rm /staging/hiseq_raw/A00977/" + x + "/analyzing.txt")
                    message = "Analyse LowPass terminée"
                    sendemail(emailto,runname,message,file)
                    #os.system("echo \"Analyse LowPass terminée\" | mutt -s \"Run Lowpass " + runname + "\" -e 'my_hdr From:LowPass Bioinfo <from@address.com>' -e 'set sendmail=/usr/sbin/ssmtp' -- luc.marchand.hsj@ssss.gouv.qc.ca")
                #os.system("cd $HOME/" + runname + " && script -c \"python /staging2/soft/CHUSJ-utils/franklin12.py " + runname + "\" -f progress.txt")
                print("Fin de la pipeline")

            elif len(glob.glob("/staging/hiseq_raw/A00977/" + x + "/LowPass*.csv")) > 1:
                message = "Il y a plus qu'une SampleSheet dans la run"
                sendemail(emailto,runname,message,file)
                #os.system("echo \"Il y a plus qu'une SampleSheet dans la run\" | mutt -s \"Run Lowpass " + x + "\" -e 'my_hdr From:LowPass Bioinfo <from@address.com>' -- luc.marchand.hsj@ssss.gouv.qc.ca")


"""
        elif len(glob.glob("/staging/hiseq_raw/LH00336/" + x + "/LowPass*.csv")) > 1:
            os.system("echo \"Il y a plus qu'une SampleSheet dans la run\" | mutt -s \"Run Lowpass " + runname + "\" -e 'my_hdr From:LowPass Bioinfo <from@address.com>' -- luc.marchand.hsj@ssss.gouv.qc.ca")
"""
