#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from datetime import datetime
from Bio import SeqIO
from Bio import Seq
import os
import sys
import yaml
import argparse
import eep_modules as mod
import eep_utils as utils
import subprocess
from itertools import product
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

#on verifie que le bon module python est chargé avant de continuer
if not mod.isloaded("python/3.6.9"):
    print("python 3.6.9 is not loaded, use: module load python/3.6.9")
    sys.exit(-1)

__version__='1.0'
args=None
config={}

#chargement des modules ici
def load_modules():
    stdout_print("Load modules:")
    moduleList=["ART/05_06_2016","repeatmasker/4.0.9_p2","bedtools/2.29","phyluce/1.6.8","samtools/1.10","stampy/1.0.32","ucsc_fatotwobit/377"]
    for module in moduleList:
        stdout_print("\t--> {}".format(module))
        mod.load(module)
        
#Fonction principale du script    
def run(ArgsVal):
    global args
    global config
    description = """ Search Ultra Conservated Elements """
    
    parser = argparse.ArgumentParser(prog="search_UCE.py.py",description=description)

    parser.add_argument("-V", "--verbose", help="full verbose", action="store_const", const="-V", default="")
    parser.add_argument("-c", "--checks", help=argparse.SUPPRESS, action="store_const", const="-c", default="")
    parser.add_argument("-t", "--test", help=argparse.SUPPRESS, action="store_const", const="-t", default="")
    parser.add_argument("-f", "--force", action="store_const", const="-f", default="")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)
    parser.add_argument("-o","--outputdir", help="output directory", required = True)
    parser.add_argument("-g","--genomes", help="genomes fasta files",nargs="*", required = True)
    parser.add_argument("-b","--basegenome", help="base genome fasta file", required = True)
    parser.add_argument("-m","--maskgenomes", help="Mask genome repetitions using repeatmasker", action="store_const", const=True, default=False)
    parser.add_argument("-s","--genspecies", help="txt file with correspondance between genome name:specie", required = True)
    parser.add_argument("--config", help="config file", default="config.yaml")
    args = parser.parse_args(ArgsVal)

    config = yaml.load(open(args.config))

    utils.setMaxParallelJobsToConfig(config['MaxParallelJobs'])
    
    datetimeNow = datetime.now()
    currentDateTime = datetimeNow.strftime("%d/%m/%Y %H:%M:%S")
    config['filesdate']=datetimeNow.strftime("%d_%m_%Y-%H_%M_%S")
    config['log_file'] = "Search_UCE_Log_{dt}.txt".format(dt=config['filesdate'])

    maskgenomes=args.maskgenomes
    genspecies=args.genspecies
    genomes=args.genomes
    basegenome=args.basegenome
    
    #paths
    config['outfolder']=args.outputdir
    config['log_folder']=os.path.join(config['outfolder'],config['logfolderName'])
    config['genomeDir']=os.path.join(config['outfolder'],'genomes')
    config['readsDir']=os.path.join(config['outfolder'],'reads')
    config['basegenomeDir']=os.path.join(config['outfolder'],'base')
    config['alignmentsDir']=os.path.join(config['outfolder'],'alignments')
    config['allAlignmentsDir']=os.path.join(config['alignmentsDir'],'all')
    config['bedDir']=os.path.join(config['outfolder'],'bed')
    config['reportsDir']=os.path.join(config['outfolder'],'reports')
    config['conservedbedfastaDir']=os.path.join(config['outfolder'],'conserved_bed_fasta')
    config['sqliteDir']=os.path.join(config['outfolder'],'sqliteDB')
    config['txtreport']=os.path.join(config['reportsDir'],'{}_report'.format(config['filesdate']))

    bedconffile=os.path.join(config['bedDir'],"bed-files.conf")
    
    if not os.path.exists(config['outfolder']):
        os.mkdir(config['outfolder'])

    stdout_print(formatTitle("Search Ultra Conservated Elements -- v{}".format(__version__)))
    stdout_print("command line: {}".format(" ".join(sys.argv)))

    #Check input files exists
    if not check_inputfiles(genomes+[genspecies],basegenome):
        stdout_print("FILES ERROR -- Exit")
        sys.exit(-1)
        
    stdout_print(formatTitle("Configurations:"))
    for param,value in config.items():
        stdout_print("\t{} : {}".format(param,value))
    #verifie que module est bien chargé
    if mod.isactivated():
        load_modules()
    else:
        stdout_print("module n'est pas chargé")

    #create all directories
    dirFile_exist(config['genomeDir'],1)
    dirFile_exist(config['readsDir'],1)
    dirFile_exist(config['basegenomeDir'],1)
    dirFile_exist(config['alignmentsDir'],1)
    dirFile_exist(config['allAlignmentsDir'],1)
    dirFile_exist(config['bedDir'],1)
    dirFile_exist(config['sqliteDir'],1)
    dirFile_exist(config['conservedbedfastaDir'],1)
    dirFile_exist(config['reportsDir'],1)

    printToReportFile(formatTitle("UCE Report {}".format(currentDateTime)))
    printToReportFile("\n")
    
    # *****************
    # * MAIN WORKFLOW *
    # *****************
    
    #Normalise fna-fasta files
    genomeFolderDic=normalise_fasta(genomes,config['genomeDir'])
    #Mask genome if paramater -m is set
    if maskgenomes:
        genomeFolderDic=mask_genomes(genomeFolderDic,genspecies)
    #Convert genomes to 2bit format
    genomeFolderDic=fastaTo2bit(genomeFolderDic)

    #reads simulation and format
    genomeFolderDic=reads_simulation(genomeFolderDic)
    genomeFolderDic=reads_format(genomeFolderDic)
    
    #base genome preparation
    basegenomeName,basegenomeFolder=prepare_base_genome(basegenome,genomeFolderDic,genspecies)
    printToReportFile("\nBase genome : {}".format(basegenome))
    
    sqldbfile=os.path.join(config['sqliteDir'],"all-to-{}.sqlite".format(os.path.basename(basegenomeFolder)))

    #simulated reads alignment against base genome
    genomeFolderDic=do_alignemnts(genomeFolderDic,basegenomeFolder)

    #Convert reads alignemnts data to coordinates in bed format
    bams_to_beds(genomeFolderDic,basegenomeFolder,bedconffile)

    #Construct Database of locus presence in genomes
    locus_presence_genomes(bedconffile,basegenomeName,sqldbfile)

    #Determining shared and conserved loci
    outconsbedfastaList=determine_shared_loci(sqldbfile,basegenomeName)

    #Extract Fasta sequences from shared an conserved loci datas
    extract_fasta_from_bedloci(outconsbedfastaList,basegenomeFolder)

    #Analyse fasta files to print stats and remove low complexity sequences
    filteredFastaList=analyse_fasta_files(outconsbedfastaList)

    merge_regionsfasta(filteredFastaList,os.path.join(config['conservedbedfastaDir'],"{}_merge_Filtered.fasta".format(basegenomeName)))
    
def merge_regionsfasta(fastalist,outputmergedfasta):
    stdout_print(formatTitle("Merge filtered fasta to {}:".format(outputmergedfasta)))
    printToReportFile("\nMerge filtered fasta to {}:".format(outputmergedfasta))
    fastadict={}
    s=0
    duplicate=0
    for n in sorted(fastalist.keys(),reverse=True):
        fastafile=fastalist[n]
        fasta=SeqIO.parse(fastafile,'fasta')
        induplicate=0
        for rec in fasta:
            seqid=str(rec.description)
            baseid=seqid.split('|')[1]
            if baseid not in fastadict.keys():
                newid="slice_{} |{}|+{}".format(s,baseid,n)
                s+=1
                fastadict[baseid]={'id':newid,'seq':str(rec.seq)}
            else:
                stdout_print("\tDuplicate sequence \"{}\" found in {}".format(baseid,fastafile))
                duplicate+=1
                induplicate+=1
        printToReportFile("\tfile:{} -- duplicates: {}".format(os.path.basename(fastafile),induplicate))
        stdout_print("\tfile:{} -- duplicates: {}".format(os.path.basename(fastafile),induplicate))
    printToReportFile("\tTotal duplicates: {}".format(duplicate))
    stdout_print("\tTotal duplicates: {}".format(duplicate))
    outfasta=""
    for seq in fastadict.values():
        outfasta+=">{}\n{}\n".format(seq['id'],seq['seq'])
    stdout_print("\tMerged fasta sequences count: {}".format(len(fastadict.keys())))
    printToReportFile("\tMerged fasta sequences count: {}".format(len(fastadict.keys())))
    with open(outputmergedfasta,'w') as outf:
        outf.write(outfasta)
        
def stats_pdf_out(outpdf,results,excludeSeq,fasta):
    datasDict={'lenDatas':[],'ShEdatas':[],'freqDatas':[],'FlenDatas':[],'FShEdatas':[],'FfreqDatas':[]}
    nseq=0
    nfilteredseq=0
    for seqid,val in results.items():
        #All datas
        datasDict['lenDatas'].append(val['len'])
        datasDict['ShEdatas'].append(val['ShE'])
        datasDict['freqDatas'].append(val['repeat']['freq'])
        nseq+=1
        #filtered datas
        if val['repeat']['freq']<config['repeat_freq_thrld']:
            datasDict['FlenDatas'].append(val['len'])
            datasDict['FShEdatas'].append(val['ShE'])
            datasDict['FfreqDatas'].append(val['repeat']['freq'])
            nfilteredseq+=1
    printToReportFile("\t{infa}:\n\t\ttotal count: {nbr}\tfiltered count: {fnbr}\ttotal length: {tlen} bp\tfiltered length: {flen} bp\tfiltered mean length: {mlen} bp\tfiltered mean entropy: {mShE}\tfiltered mean base repetition freq: {mfreq}".format(infa=os.path.basename(fasta),nbr=nseq,fnbr=nfilteredseq,tlen=sum(datasDict['lenDatas']),flen=sum(datasDict['FlenDatas']),mlen=round(np.mean(datasDict['FlenDatas']),2),mShE=round(np.mean(datasDict['FShEdatas']),2),mfreq=round(np.mean(datasDict['FfreqDatas']),2)))
    titles={'lenDatas':'Sequences Length distribution (n={})'.format(nseq),'ShEdatas':'Sequences entropy distribution','freqDatas':'max motifs frequence distribution','FlenDatas':'Filtered sequences Length distribution  (n={})'.format(nfilteredseq),'FShEdatas':'Filtered sequences entropy distribution','FfreqDatas':'Filtered max motifs frequence distribution'}
    with PdfPages(outpdf) as pdf:
        for header,datas in datasDict.items():
            title="{}\n{}".format(titles[header],os.path.basename(fasta[:fasta.rfind('.')]))
            dist=plot_histogram(datas,title)
            if header=='freqDatas':
                plt.axvline(config['repeat_freq_thrld'],0,max(dist[0]),label='repeat frequence threshold\nfor filtered fasta',color='red')
                plt.legend()
            pdf.savefig(plt.gcf())
    
def plot_histogram(datas,title):
    fig=plt.figure()
    plt.clf()
    dist=plt.hist(datas,bins=50)
    plt.title(title)
    return dist
    
def analyse_fasta_files(outconsbedfastaList):
    stdout_print(formatTitle("Analyse Fasta files:"))
    printToReportFile("\nAnalyse Fasta files:")
    filteredFastaList={}
    for bed,fasta,n in outconsbedfastaList:
        filteredFastaName="{}_Filtered.fasta".format(fasta[:fasta.rfind('.')])
        results,excludeSeq=find_repeated_sequence(fasta,filteredFastaName)
        filteredFastaList[n]=filteredFastaName
        stdout_print("\n\t=> fasta file {}".format(fasta))
        stdout_print("\t\tSeq Id\tlength\tShE\tcodon\tFreq\tsequence")
        for seqid,val in excludeSeq.items():
            stdout_print("\t\t{sid}\t{seqlen}\t{she}\t{codonfreq}\t{freqmax}\t{tseq}...".format(sid=seqid,seqlen=val['len'],she=val['ShE'],codonfreq=val['codonfreq'],freqmax=val['freqmax'],tseq=val['tseq']))
        outpdf=os.path.join(config['reportsDir'],"{dt}_distributions_{fas}.pdf".format(fas=os.path.basename(fasta[:fasta.rfind('.')]),dt=config['filesdate']))
        stats_pdf_out(outpdf,results,excludeSeq,fasta)
    return filteredFastaList
        
def find_repeated_sequence(fastafile,outnewfasta=None):
    rRange=config['repeat_range']
    fasta = SeqIO.parse(fastafile,'fasta')
    seqstats={}
    newfasta=""
    excludefasta=""
    excludeSeq={}
    for rec in fasta:
        sequence=str(rec.seq).upper()
        seqstats[str(rec.id)]={}
        seqstats[str(rec.id)]['len']=len(sequence)
        seqstats[str(rec.id)]['ShE']=get_ShannonEntropy(sequence)
        seqstats[str(rec.id)]['repeat']={}
        freqmax=0
        rslt_fremax_codon={}
        for r in range(rRange[0],rRange[0]+1):
            repet_results=dict(("".join(codon), sequence.upper().count(''.join(codon))) for codon in product("ATCG", repeat=r))
            for codon,nbr in repet_results.items():
                freq=float(nbr)/(float(len(sequence))/r)
                if freq>freqmax:
                    freqmax=freq
                    rslt_fremax_codon={'codon':codon,'n':nbr,'freq':freq}
        seqstats[str(rec.id)]['repeat']=rslt_fremax_codon
        if freqmax<config['repeat_freq_thrld']:
            newfasta+=">{}\n{}\n".format(rec.description,sequence)
        else:
            excludeSeq[str(rec.id)]={'len':len(sequence),'tseq':sequence[:30],'ShE':seqstats[str(rec.id)]['ShE'],'freqmax':rslt_fremax_codon['freq'],'codonfreq':rslt_fremax_codon['codon']}
            excludefasta+=">{}\n{}\n".format(rec.description,sequence)
    if outnewfasta:
        with open(outnewfasta,'w') as nwf:
            nwf.write(newfasta)
        with open("{}.excluded".format(outnewfasta),'w') as exf:
            exf.write(excludefasta)
    return seqstats,excludeSeq

def get_ShannonEntropy(dna):
    ShEntropy = 0
    Udna = dna.upper()
    alphabet=get_ambiguous_dna("N")

    for base in alphabet:
        freq = float(Udna.count(base))/float(len(Udna))
        if freq==0: freq = 0.00000001
        ShEntropy = ShEntropy + freq*math.log(freq,2)
    return -ShEntropy

def get_ambiguous_dna(seq):
    d = Seq.IUPAC.IUPACData.ambiguous_dna_values
    return list(map("".join, product(*map(d.get, seq))))

def extract_fasta_from_bedloci(outconsbedfastaList,basegenomeFolder):
    basegenome2bit="{}.fasta.2bit".format(basegenomeFolder)
    for rec in outconsbedfastaList:
        cmd=config['extractFasta_cmd'].format(bedfile=rec[0],gbase2bit=basegenome2bit,outfasta=rec[1])
        job = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = job.communicate()
        stdout_print("\t{}".format(cmd))
        stdout_print("\t\t=> {}".format(stdout.decode("utf-8")))
        
def determine_shared_loci(sqldbfile,basegenomeName):
    outshared={}
    outconsbedfastaList=[]
    stdout_print(formatTitle("Determining shared, conserved, loci:"))
    printToReportFile("\nDetermining shared, conserved, loci:")
    cmd="phyluce_probe_query_multi_merge_table --db {} --base-taxon {}".format(sqldbfile,basegenomeName)
    job = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = job.communicate()
    stdout_print("\t{}".format(cmd))
    for line in stdout.decode("utf-8").split('\n'):
        stdout_print("\t\t{}".format(line))
        printToReportFile("\t{}".format(line))
        if "+" in line:
            tmp=line.split("+")[1].strip().split()
            outshared[int(tmp[0])]=float(tmp[2].replace(',',''))
    for n in range(len(outshared)-config['number_most_conserved_loci'],len(outshared)):
        outconsbed=os.path.join(config['conservedbedfastaDir'],"{gbase}+{nbr}.bed".format(gbase=basegenomeName,nbr=n))
        outconsfasta=os.path.join(config['conservedbedfastaDir'],"{gbase}+{nbr}.fasta".format(gbase=basegenomeName,nbr=n))
        outconsbedfastaList.append((outconsbed,outconsfasta,n))
        cmd="phyluce_probe_query_multi_merge_table --db {sqldb} --base-taxon {gbase} --output {bedfile} --specific-counts {nbr}".format(sqldb=sqldbfile,gbase=basegenomeName,bedfile=outconsbed,nbr=n)
        stdout_print("\t{}".format(cmd))
        printToReportFile("\tcounts for \"{gbase}+{nbr}\"".format(gbase=basegenomeName,nbr=n))
        job = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = job.communicate()
        dictStr=stdout.decode("utf-8")
        dictStr=eval(dictStr[dictStr.find('{'):dictStr.rfind(')')])
        for gen,nbr in dictStr.items():
            stdout_print("\t\t{}\t{}".format(nbr,gen))
            printToReportFile("\t\t{}\t{}".format(nbr,gen))
    return outconsbedfastaList

def locus_presence_genomes(bedfileconf,basegenomeName,sqldbfile):
    stdout_print(formatTitle("Determining locus presence in multiple genomes:"))
    if not os.path.exists(sqldbfile):
        cmd="phyluce_probe_get_multi_merge_table --conf {bedconf} --base-taxon {gbase} --output {sqldb}".format(bedconf=bedfileconf,gbase=basegenomeName,sqldb=sqldbfile)
        stdout_print("\t{}".format(cmd))
        os.system(cmd)
    else:
        stdout_print("\t{} already exists".format(sqldbfile))
        
def bams_to_beds(ingenomeDic,basegenomeFolder,bedconffile):
    stdout_print(formatTitle("Bams to beds convertion:"))
    printToReportFile("Bed informations:")
    bedfileconf="[beds]"
    for gname,values in ingenomeDic.items():
        bedfile="{}.bed".format(os.path.join(config['bedDir'],os.path.basename(values['allbam'])))
        sortbedfile="{}.sort.bed".format(os.path.join(config['bedDir'],os.path.basename(values['allbam'])))
        mergedbedfile="{}.merge.bed".format(os.path.join(config['bedDir'],os.path.basename(values['allbam'])))
        stripbedfile="{}.strip.bed".format(os.path.join(config['bedDir'],os.path.basename(values['allbam'])))
        bedfileconf+="\n{}:{}".format(values['genomeName'],stripbedfile)
        printToReportFile("\n\tfor genome {}".format(gname))
        if not os.path.exists(stripbedfile):
            cmd="bedtools bamtobed -i {inbam} -bed12 > {bedname}".format(inbam=values['allbam'],bedname=bedfile)
            os.system(cmd)
            stdout_print("\t{} => {}".format(gname,cmd))
            cmd="bedtools sort -i {inbed} > {sortbed}".format(inbed=bedfile,sortbed=sortbedfile)
            os.system(cmd)
            stdout_print("\t{} => {}".format(gname,cmd))
            cmd="bedtools merge -i {insortbed} > {mergebed}".format(insortbed=sortbedfile,mergebed=mergedbedfile)
            os.system(cmd)
            stdout_print("\t{} => {}".format(gname,cmd))
            cmd=config['removeRepetInterval_cmd'].format(mergedbed=mergedbedfile,basegenome2bit="{}.fasta.2bit".format(basegenomeFolder),outstripbed=stripbedfile)
            job = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = job.communicate()
            stdout_print("\tphyluce_probe_strip_masked_loci_from_set out => {}".format(stdout.strip().decode("utf-8")))
            stdout_print("\t{} => {}".format(gname,cmd))
        else:
            stdout_print("\tBeds files already exists for {}".format(gname))
            
        report_bedFile(bedfile)
        report_bedFile(mergedbedfile)
        report_bedFile(stripbedfile)
        
        cmd="wc -l {}".format(mergedbedfile)
        job = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = job.communicate()
        stdout_print("\t\t=> {} contains {} merged regions".format(os.path.basename(mergedbedfile),int(stdout.strip().split()[0])))
        cmd="wc -l {}".format(stripbedfile)
        job = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = job.communicate()
        stdout_print("\t\t=> {} contains {} regions".format(os.path.basename(stripbedfile),int(stdout.strip().split()[0])))
        
    with open(bedconffile,'w') as bf:
        bf.write(bedfileconf)

def report_bedFile(bedfile):
    regionsSize=[]
    chrname=set()
    with open(bedfile,'r') as bfile:
        for line in bfile:
            tmp=line.strip().split()
            chrv=tmp[0]
            pstart=int(tmp[1])
            pend=int(tmp[2])
            regionsSize.append(abs(pend-pstart))
            chrname.add(chrv)
    printToReportFile("\t\t=> {bedname}\tChr number: {chrnbr}\tregion min size: {mins} bp\tregion max size:{maxs} bp\tregion mean size: {means} bp".format(bedname=os.path.basename(bedfile),chrnbr=len(chrname),mins=min(regionsSize),maxs=max(regionsSize),means=round(np.mean(regionsSize),2)))
    
def do_alignemnts(ingenomeDic,basegenomeFolder):
    stdout_print(formatTitle("Align reads to the base genome:"))
    for gname,values in ingenomeDic.items():
        outbam=os.path.join(config['alignmentsDir'],"{}-to-{}-MAPPING.bam".format(gname,os.path.basename(basegenomeFolder)))
        ingenomeDic[gname]['bam']=outbam
        ingenomeDic[gname]['allbam']=os.path.join(config['allAlignmentsDir'],os.path.basename(outbam))
        if not os.path.exists(outbam):
            cmd=config['align_cmd'].format(gbase=basegenomeFolder,read=values['read'],outbam=outbam)
            stdout_print("\t{}".format(cmd))
            os.system(cmd)
            stdout_print("\tcreate link FROM \"{}\" TO \"{}\"".format("../{}".format(os.path.basename(outbam)),os.path.join(config['allAlignmentsDir'],os.path.basename(outbam))))
            os.symlink("../{}".format(os.path.basename(outbam)),os.path.join(config['allAlignmentsDir'],os.path.basename(outbam)))
        else:
            stdout_print("\t{} already exists".format(outbam))
    return ingenomeDic
        
def check_inputfiles(infiles,basegenome):
    stdout_print(formatTitle("Check input files:"))
    check=True
    for f in infiles:
        if not os.path.exists(f):
            stdout_print("\tERROR: {} not found".format(f))
            check=False
    if basegenome not in infiles:
        stdout_print("\tERROR: basegenome {} not in genome input list".format(basegenome))
        check=False
        
    return check

def prepare_base_genome(basegenome,genomeFolderDic,genomespeciesfile):
    stdout_print(formatTitle("Prepare base genome:"))
    gspec=get_genomes_species(genomespeciesfile)
    basegenomeFasta=os.path.basename(basegenome)
    basegenomeFolder=basegenomeFasta[:basegenomeFasta.rfind('.')]
    outbasegenome=os.path.join(config['basegenomeDir'],basegenomeFasta)
    outbasegenome2bit=os.path.join(config['basegenomeDir'],"{}.2bit".format(basegenomeFasta))
    basegenomeFastaFromDic=genomeFolderDic[basegenomeFolder]['fasta']
    stdout_print("\tcreate link FROM \"{}\" TO \"{}\"".format(basegenomeFastaFromDic[basegenomeFastaFromDic.find('/'):],outbasegenome))
    if not os.path.islink(outbasegenome):
        os.symlink("..{}".format(basegenomeFastaFromDic[basegenomeFastaFromDic.find('/'):]),outbasegenome)
    if not os.path.islink(outbasegenome2bit):
        os.symlink("..{}.2bit".format(basegenomeFastaFromDic[basegenomeFastaFromDic.find('/'):]),outbasegenome2bit)
    cmd="stampy.py --species=\"{specie}\" --assembly=\"{asm}\" -G {asm} {fasta}".format(specie=gspec[basegenomeFolder],asm=os.path.join(config['basegenomeDir'],basegenomeFolder),fasta=outbasegenome)
    stdout_print("\t{}".format(cmd))
    os.system(cmd)
    cmd="stampy.py -g {asm} -H {asm}".format(asm=os.path.join(config['basegenomeDir'],basegenomeFolder))
    stdout_print("\t{}".format(cmd))
    os.system(cmd)
    return basegenomeFolder.replace(".","_"),os.path.join(config['basegenomeDir'],basegenomeFolder)
    
def reads_format(ingenomeDic):
    stdout_print(formatTitle("Reads concatenates and compression:"))
    jblst = utils.Jobslist("gzip")
    for gname,values in ingenomeDic.items():
        reads=values['readsfw']
        outconcatreads=os.path.join(config['readsDir'],"{}-reads.fq".format(gname))
        outconcatreadsgz="{}.gz".format(outconcatreads)
        ingenomeDic[gname]['read']=outconcatreadsgz
        if not os.path.exists(outconcatreadsgz):
            stdout_print("\t--> concat reads : {}".format(reads))
            with open(outconcatreads,'w') as outconcat:
                for read in reads:
                    with open(read,'r') as inrd:
                        outconcat.write(inrd.read().strip())
                        outconcat.write("\n")
                    os.remove(read)
            cmd="gzip {}".format(outconcatreads)
            stdout_print("\t--> compress: {}".format(cmd))
            jblst.add_a_job(cmd,cmd,outconcatreadsgz)
        else:
            stdout_print("\t--> reads are already concatenated and compressed: {}".format(gname))
    utils.trun(args, jblst)
    return ingenomeDic
        
def reads_simulation(ingenomeDic):
    jblst = utils.Jobslist("artillumina")
    stdout_print(formatTitle("Reads simulation:"))
    for gname,values in ingenomeDic.items():
        readsprefix=os.path.join(config['readsDir'],"{}-reads".format(gname))
        outconcatreads=os.path.join(config['readsDir'],"{}-reads.fq".format(gname))
        outconcatreadsgz="{}.gz".format(outconcatreads)
        cmd=config['artillumina_cmd'].format(infasta=values['fasta'],readsprefix=readsprefix)
        target="{}1.fq".format(readsprefix)
        ingenomeDic[gname]['readsfw']=("{}1.fq".format(readsprefix),"{}2.fq".format(readsprefix))
        if not os.path.exists(target) and not os.path.exists(outconcatreadsgz):
            stdout_print("\t--> artillumina simulation: {}".format(cmd))
            jblst.add_a_job(cmd,"artillumina {}".format(gname),target)
        else:
            stdout_print("\t--> {} READS Already simulated".format(gname))
    utils.trun(args, jblst)
    return ingenomeDic
        
def fastaTo2bit(ingenomeDic):
    stdout_print(formatTitle("Convert genomes to 2bit format:"))
    for gname,values in ingenomeDic.items():
        fastaFile=values['fasta']
        out2bitFile="{}.2bit".format(fastaFile)
        ingenomeDic[gname]['2bit']=out2bitFile
        if not os.path.exists(out2bitFile):
            cmd="faToTwoBit {infasta} {out2bit}".format(infasta=fastaFile,out2bit=out2bitFile)
            stdout_print("\t--> 2bit convert: {}".format(cmd))
            os.system(cmd)
        else:
            stdout_print("\t--> \"{}\" Aldready 2bit converted".format(fastaFile))
    return ingenomeDic
            
def normalise_fasta(fnalist,outputgenomes):
    genomeFolderDic={}
    stdout_print(formatTitle("Normalize fasta genomes:"))
    printToReportFile("Genomes informations:")
    for fna in fnalist:
        fname=os.path.basename(fna)
        fname=fname[:fname.rfind(".")]
        genomefolder=os.path.join(outputgenomes,fname)
        genomeName=os.path.basename(genomefolder).replace(".","_")
        dirFile_exist(genomefolder,1,False)
        outfasta=os.path.join(genomefolder,"{}.fasta".format(fname))
        genomeFolderDic[fname]={'folder':genomefolder,'fasta':outfasta,'genomeName':genomeName}
        genomesize=0
        sequenceNbr=0
        maskedBp=0
        if not os.path.exists(outfasta):
            stdout_print("\t--> START Normalise {}".format(fna))
            with open(fna, "rU") as infile:
                with open(outfasta, "w") as outf:
                    for seq in SeqIO.parse(infile, 'fasta'):
                        seq.name = ""
                        seq.description = ""
                        genomesize+=len(seq.seq)
                        maskedBp+=str(seq.seq).upper().count("N")
                        sequenceNbr+=1
                        outf.write(seq.format('fasta'))
        else:
            for seq in SeqIO.parse(fna, 'fasta'):
                genomesize+=len(seq.seq)
                sequenceNbr+=1
                maskedBp+=str(seq.seq).upper().count("N")
            stdout_print("\t--> {} ALREADY NORMALIZED".format(fna))
        printToReportFile("\n\t=> {gen} -- sequences number: {snb} -- size: {size} bp -- N bases: {msize} bp ({pmsize} %)".format(gen=os.path.basename(fna),snb=sequenceNbr,size=genomesize,msize=maskedBp,pmsize=round(float(maskedBp)/float(genomesize),2)))
    return genomeFolderDic

def get_genomes_species(genomespeciesfile):
    gspec={}
    with open(genomespeciesfile,'r') as gspecfile:
        for line in gspecfile:
            tmp=line.strip().split(',')
            if len(tmp)==2:
                gspec[tmp[0]]=tmp[1]
    return gspec

def mask_genomes(genomeFolderDic,maskgenomes):
    results={}
    stdout_print(formatTitle("Mask genome repetitions:"))
    maskInfos=get_genomes_species(maskgenomes)
    for gname,values in genomeFolderDic.items():
        maskedFastaFile="{}.masked".format(values['fasta'])
        results[gname]={'folder':values['folder'],'fasta':maskedFastaFile,'genomeName':values['genomeName']}
        if not os.path.exists(maskedFastaFile):
            cmd="RepeatMasker -pa {ncore} -dir {outdir} -species \"{specie}\" {fasta} > {rpmasker_out}".format(ncore=config['repeatmasker_cores'],outdir=values['folder'],specie=maskInfos[gname],fasta=values['fasta'],rpmasker_out="{}/repeatMasker_output".format(values['folder']))
            stdout_print("\t--> RepeatMasker: {}".format(cmd))
            os.system(cmd)
        else:
            stdout_print("\t--> RepeatMasker already done for {}".format(gname))
    return results

def printToReportFile(txt):
    with open(config['txtreport'],'a+') as reportFile:
        reportFile.write("{}\n".format(txt))

def formatTitle(txt):
    starline="*"*(len(txt)+4)
    outitle="{star}\n* {title} *\n{star}\n".format(title=txt,star=starline)
    return outitle

def addTextToLogFile(logtxt,printTime=True):
    if printTime:
        ttime = "{}\t".format(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]"))
    else:
        ttime=""
        
    logFolder = config['log_folder']
    
    dirFile_exist(logFolder,1,False)

    logFile = os.path.join(logFolder,config['log_file'])
    lf = open(logFile,'a+')
    if logtxt.count('\n')>0:
        for l in logtxt.split('\n'):
            lf.write("{tm}{log}\n".format(tm=" "*len(ttime),log=l))
    else:
        lf.write("{tm}{log}\n".format(tm=ttime,log=logtxt))
    lf.close()

def stdout_print(txt,printInLog=True):
    print(txt)
    if printInLog:
        addTextToLogFile(txt)
                
#Test if a file or directory exist if creatdir set to 0: exit
# if creatdir set to 1: create new directory
def dirFile_exist(path,creatdir = 0,printInLog=True):
    if not os.path.exists(path):
        if (creatdir == 0):
            exitProg("ERROR: {} not found".format(path))
        else:
            dirmsg = "-- create {} directory".format(path)
            stdout_print(dirmsg,printInLog)
            os.makedirs(path)

def exitProg(msg=''):
    if msg!='':
        addTextToLogFile("PROGRAM EXIT\t{}".format(msg))
    sys.exit(msg)
        
if __name__=="__main__":
    ArgsVal = sys.argv[1:]
    run(ArgsVal)
