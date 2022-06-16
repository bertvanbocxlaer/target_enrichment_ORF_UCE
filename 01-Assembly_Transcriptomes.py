#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import utils
import time
import threading
import argparse
import sys
import os
from datetime import datetime
import subprocess
import eep_modules as mod
#on verifie que le bon module python est chargé avant de continuer
if not mod.isloaded("python/2.7.15"):
    print("python 2.7.15 is not loaded, use: module load python/2.7.15")
    sys.exit(-1)
    
__version__="0.0.2"

args = None
config = {}

def run(ArgqsVal):
    
    load_modules()
    global args
    global config
    description = """ Agalma python pipeline Assembly """

    create_dir("./agalma_DB")
    
    parser = argparse.ArgumentParser(prog="Assemble_Transcriptomes.py",description=description)

    parser.add_argument('-ver', '--version', action='version', version='%(prog)s ' + __version__)
    parser.add_argument("-V", "--verbose", help="full verbose", action="store_const", const="-V", default="")
    parser.add_argument("-v", "--tinyverbose", help="verbose",action="store_const", const="-v", default="")
    parser.add_argument("-t", "--test", help=argparse.SUPPRESS, action="store_const", const="-t", default="")
    parser.add_argument("-f", "--force", action="store_const", const="-f", default="")
    parser.add_argument("-c", "--checks", help=argparse.SUPPRESS, action="store_const", const="-c", default="")
    parser.add_argument("-q", "--doqc",help="do fastq quality control",action="store_const", const=True, default=False)
    parser.add_argument("-o","--outputdir", help="output directory", required = True)
    parser.add_argument("-r","--reads_list",help="text file with reads list",required = True)
    parser.add_argument("-i","--itis_id",help="itis ID",default=None)
    parser.add_argument("-d","--agalmaDBFld", help="Agalma Database file", default="/srv/titan/baobab/projects1/Evolink/pipeline_Assembly/agalma_DB/agalma{}.sqlite".format(get_currentTime(True)))
    args = parser.parse_args(ArgsVal)
    config['outputdir'] = args.outputdir
    create_dir(config['outputdir'])

    config['logfilename'] = os.path.join(config['outputdir'],"log_pipeline_Assembly_{}.txt".format(get_currentTime(True)))
    add_Text_ToLog("Assembly Pipeline with Agalma Start")

    export_AgalmaDB(args.agalmaDBFld)

    samplesDic = parse_readsFile(args.reads_list)

    nbNotFound = check_fqFilesExists(samplesDic)
    if nbNotFound>0 :
        add_Text_ToLog("EXIT: {} FILE(S) NOT FOUND".format(nbNotFound))
        sys.exit(2)
    
    #Workflow for Agalma pipeline
    insert_AgalmaCatalog(samplesDic,args.itis_id)
    if args.doqc:
        fastq_QC(samplesDic)
    do_assembly(samplesDic)
    do_reports(samplesDic)
    add_Text_ToLog("Assembly Pipeline with Agalma ENDED...")
    send_a_mail("transcriptome assembly ENDED {}".format(get_currentTime()),"Assemblages Claudia","mathieu.genete@univ-lille.fr")

#chargement des modules ici
def load_modules():
        mod.load("agalma/2.0.0")
        
def send_a_mail(message,sujet,dest):
    cmd = 'echo "{}" | mail -s "{}" {}'.format(message,sujet,dest)
    os.system(cmd)
    
def check_fqFilesExists(samplesDic):
    add_Text_ToLog("Check fastq files path:")
    nbNotFound=0
    for ind,values in samplesDic.items():
        for fq in [values['fwd'],values['rev']]:
            if not os.path.exists(fq):
                add_Text_ToLog("\t=> {} NOT FOUND for sample {}".format(fq,ind))
                nbNotFound+=1
    if nbNotFound==0:
        add_Text_ToLog("\tAll fq files are OK")
    return nbNotFound  
                
            
def progress_assembly(samplesDic,logfname):
    nbrSample = len(samplesDic)
    indivTerminated=set()
    inProgress=set()
    print("--- ok ---")
    add_Text_ToLog("Launch progress assembly function for {} samples".format(nbrSample),True,logfname)
    while len(indivTerminated)<nbrSample:
        time.sleep(2)
        for ind,values in samplesDic.items():
            if test_folderFilePresent("transcriptome-",values['outfld']) and ind not in inProgress:
                inProgress.add(ind)
                add_Text_ToLog("PROGRESS ASSEMBLY: assembly started for {}... {}/{}".format(ind,len(inProgress),nbrSample),True,logfname)
            if test_file_inFolderSubName("transcriptome-",values['outfld'],"{}.fa".format(ind)) and ind not in indivTerminated:
                indivTerminated.add(ind)
                add_Text_ToLog("PROGRESS ASSEMBLY: assembly ended for {}... {}/{}".format(ind,len(indivTerminated),nbrSample),True,logfname)
    add_Text_ToLog("STOP Launch progress assembly function...".format(nbrSample),True,logfname)
        
def fastq_QC(samplesDic):
    add_Text_ToLog("Quality control fastq:")
    jblst = utils.Jobslist("Agalma qc")
    tmpdir = os.getcwd()
    for ind,values in samplesDic.items():
        cmd="cd {}; agalma qc -i {}".format(values['outfld'],ind)
        if not test_folderFilePresent("qc-",values['outfld']):
            add_Text_ToLog("\t=> do QC for {} : {}".format(ind,cmd))
            jblst.add_a_job(cmd,"Agalma qc","")
        else:
            add_Text_ToLog("\t{} qc already present. skip...".format(ind))

    if jblst.get_SizeOfJoblist()>0:
        utils.trun(args, jblst, False)
    os.chdir(tmpdir)

def do_assembly(samplesDic):
    utils.setMaxParallelJobsToConfig(1)
    add_Text_ToLog("Launch Assembly:")
    jblst = utils.Jobslist("Agalma qc")
    tmpdir = os.getcwd()
    for ind,values in samplesDic.items():
        #modif du 21/01 pour analyse de Mol04 et Mol08
        #cmd="cd {}; agalma transcriptome --skip setup_rrna,filter_data -i {}".format(values['outfld'],ind)
        cmd="cd {}; agalma transcriptome -i {}".format(values['outfld'],ind)
        if not test_file_inFolderSubName("transcriptome-",values['outfld'],"{}.fa".format(ind)):
            jblst.add_a_job(cmd,"Agalma transcriptome","")
            add_Text_ToLog("\t=> do transcriptome for {} : {}".format(ind,cmd))
        else:
            add_Text_ToLog("\t{} transcriptome already present. skip...".format(ind))
    if jblst.get_SizeOfJoblist()>0:
        t1 = threading.Thread(target=progress_assembly, args=(samplesDic,config['logfilename'],))
        t1.start()
        utils.trun(args, jblst)
    os.chdir(tmpdir)
    add_Text_ToLog("Assembly ENDED...")

def do_reports(samplesDic):
    utils.setMaxParallelJobsToConfig(10)
    add_Text_ToLog("Launch Reports:")
    jblst = utils.Jobslist("Agalma report")
    for ind,values in samplesDic.items():
        report_path = os.path.join(values['outfld'],"REPORT_{}".format(ind))
        cmd="agalma report -i {} -o {}".format(ind,report_path)
        if not os.path.exists(report_path):
            jblst.add_a_job(cmd,"Agalma report","")
            add_Text_ToLog("\t=> do report for {} : {}".format(ind,cmd))
        else:
            add_Text_ToLog("\t{} report already present. skip...".format(ind))
    if jblst.get_SizeOfJoblist()>0:
        utils.trun(args, jblst)

    
def test_file_inFolderSubName(fldname,testpath,fileName):
    fldlist = os.listdir(testpath)
    for fld in fldlist:
        if fldname in fld:
            ispath = os.path.join(os.path.join(testpath,fld),fileName)
            if os.path.exists(ispath):
                return True
    return False

def test_folderFilePresent(fldname,testpath):
    fldlist = os.listdir(testpath)
    for fld in fldlist:
        if fldname in fld:
            return True
    return False

def add_Text_ToLog(txtline,showmsg=True,logfname=None):
    if not logfname:
        logfname=config['logfilename']
    lf = open(logfname,'a')
    logsmsg="{}\t{}".format(get_currentTime(),txtline)
    if showmsg:
        print(logsmsg)
    lf.write("{}\n".format(logsmsg))
    lf.close
    
def get_currentTime(forFile=False):
    if forFile:
        strTime = "%d_%m_%Y-%Hh%Mm%S"
    else:
        strTime = "%d/%m/%Y %H:%M:%S"
    return datetime.now().strftime(strTime)

def insert_AgalmaCatalog(samplesDic,itisID):
    if itisID!=None:
        useItisID=True
    else:
        useItisID=False

    utils.setMaxParallelJobsToConfig(10)   
    add_Text_ToLog("Insert into Agalma catalog:")
    for indiv,values in samplesDic.items():
        if useItisID:
            iID = itisID
        else:
            iID = values['itisID']
            
        if str(iID)=="0":
            cmd = "agalma catalog insert --paths {fwd} {rev} -i {sname}".format(fwd=values['fwd'],rev=values['rev'],sname=indiv)
        else:
            cmd = "agalma catalog insert --paths {fwd} {rev} --itis_id {itisID} -i {sname}".format(fwd=values['fwd'],rev=values['rev'],itisID=iID,sname=indiv)
        add_Text_ToLog("\t{}".format(cmd))
        os.system(cmd)

def create_dir(dirname):
    if not os.path.exists(dirname):
        os.mkdir(dirname)
        
def parse_readsFile(reads_list_path):
    rf = open(reads_list_path,'r')
    samplesDic={}
    l=0
    for line in rf:
        l+=1
        tmp = line.strip().split(",")
        if len(tmp)==4:
            outfld=os.path.join(config['outputdir'],tmp[0])
            samplesDic[tmp[0]]={'fwd':tmp[1],'rev':tmp[2],'outfld':outfld,'itisID':tmp[3]}
            add_Text_ToLog("Create folder {}".format(outfld))
            create_dir(outfld)
        else:
            add_Text_ToLog("ERROR READS LINE {} : {}".format(l,line.strip()))
    rf.close
    add_Text_ToLog("Number of samples: {}".format(len(samplesDic)))
    return samplesDic

def export_AgalmaDB(agalmaDB_path):
    add_Text_ToLog("export AGALMA_DB env variable: {}".format(agalmaDB_path))
    os.putenv("AGALMA_DB",agalmaDB_path)

if __name__=="__main__":
    ArgsVal = sys.argv[1:]
    run(ArgsVal)
