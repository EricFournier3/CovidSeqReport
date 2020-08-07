# coding=utf-8

"""
Eric Fournier
2020-08-06
"""

#TODO
"""
- ajouter les metrics dans la MySQL
- ne pas rouler ce script sur beluga => un script pour importer seulement les nom de plaque et csv
=> script bash ou modifier script sandrine pour creer seulement les metric.csv et ensuite on importe, ou scp on only specific extension
"""


from Bio import SeqIO
import os
import sys
import stat
from datetime import date
from Bio import SeqIO
import re
import logging
import pandas as pd
import re
from datetime import datetime
import glob



_DEBUG = True
logging.basicConfig(level=logging.DEBUG)
pd.options.display.max_columns = 100


class PlateDirManager:
    def __init__(self):
        pass
 
    @staticmethod
    def GetBaseDir(_debug):
        if _debug:
            return "/data/Users/Eric/Covid19/TestRapportSequence/PLATE_DIR"
        else:
            return ""

class PlateManager():
    def __init__(self):
        self.plates_list = []
        self.samples_list = []
    
    def AddPlate(self,plate_name):
        plate = Plate(plate_name)
        self.plates_list.append(plate)

    def GetPlates(self):
        return self.plates_list

    def GetPlatesNames(self):
        names = []
        for plate in self.plates_list:
            names.append(plate.GetName())

        return names

    def SetSampleList(self):
        for plate in self.plates_list:
            for sample in plate.GetSamplesList():
                self.samples_list.append(sample)

    def GetSamples(self):
        return(self.samples_list)
                
class Plate:
    def __init__(self,plate_name):
        self.plate_name = plate_name
        self.name = plate_name
        self.samples_manager = SamplesManager(self)

        self.SetSampleList()

    def GetName(self):
        return self.name

    def SetSampleList(self):
        self.samples_manager.SetSamplesList()

    def GetSamplesList(self):
        return(self.samples_manager.GetSamplesList())

    def GetSamplesNames(self):
        return(self.samples_manager.GetSamplesNames())

class SamplesManager:
    def __init__(self,plate):
        self.parent_plate = plate 
        self.samples_names = []
        self.samples_list = []
    
    def SetSamplesList(self):
        for spec_dir in os.listdir(os.path.join(basedir,self.parent_plate.GetName())):
            spec_name,tech,date = Utils.ParseSampleDir(spec_dir)
            if spec_name not in self.samples_names:
                self.samples_names.append(spec_name)
                new_sample = Sample(spec_name,self.GetParentPlate())
                self.samples_list.append(new_sample) 
                new_sample.UpdateSeqInfo(tech,date,spec_dir)
            else:
               sample = Utils.GetObjectByName(spec_name,self.samples_list)
               sample.UpdateSeqInfo(tech,date,spec_dir)

    def GetSamplesList(self):
        return self.samples_list

    def GetSamplesNames(self):
        return self.samples_names

    def GetParentPlate(self):
        return self.parent_plate
    
class Sample:
    def __init__(self,name,parent_plate):
        self.sample_name = name
        self.parent_plate = parent_plate

        self.metric_manager = MetricsManager(self)

    def UpdateSeqInfo(self,tech,date,spec_dir):
        date = datetime.strptime(date,'%Y%m%d') 

        self.AddNewMetric(spec_dir,date,tech)


    def AddNewMetric(self,spec_dir,date,tech):
        full_metric_path = os.path.join(PlateDirManager.GetBaseDir(_DEBUG),self.parent_plate.GetName(),spec_dir)
        
        try:
            metric_file = glob.glob(full_metric_path + '/*metrics.csv')[0]
            self.metric_manager.AddNewMetric(metric_file,date,tech) 
        except Exception as e:
            print(sys.exc_info())
            logging.warning("Aucun metric dans " + full_metric_path)

    def GetSampleName(self):
        return self.sample_name

    def GetTechAndDate(self):
        return self.tech_date

    def BuildSeqReport(self):
        for metric in self.metric_manager.GetMetrics():
            MarkdownWriter.BuildSeqReport(metric)


class MetricsManager:
    def __init__(self,sample):
        self.sample = sample
        self.metrics_list = []

        nanopore_col_renamed = {'sample':'SAMPLE','cons.perc.N':'PERC_N','cons.len':'CONS_LEN','cons.perc.GC':'PERC_GC',
                                'fq.size.pass':'FASTQ_SIZE_PASS','bam.perc.align':'BAM_PERC_ALIGN','bam.mean.cov':'BAM_MEAN_COV',
                                'bam.med.cov':'BAM_MED_COV','bam.max.min.ratio':'BAM_MAX_MIN_RATIO','bam.perc.50x':'BAM_PERC_50X',
                                'bam.perc.100x':'BAM_PERC_100X','bam.perc.250x':'BAM_PERC_250X','bam.perc.500x':'BAM_PERC_500X',
                                'bam.perc.1000x':'BAM_PERC_1000X','bam.perc.2000x':'BAM_PERC_2000X'}

        illumina_col_renamed = {'sample':'SAMPLE','cons.per.N':'PERC_N','cons.len':'CONS_LEN','cons.perc.GC':'PERC_GC',
                                'fq.trim.pass':'FASTQ_SIZE_PASS','bam.primertrim.pass':'BAM_PRIMERTRIM_PASS','bam.perc.align':'BAM_PERC_ALIGN',
                                'bam.mean.cov':'BAM_MEAN_COV','bam.med.cov':'BAM_MED_COV','bam.max-min/mean.cov':'BAM_MAXMIN_MEANCOV',
                                'bam.perc.50x':'BAM_PERC_50X','bam.perc.100x':'BAM_PERC_100X','bam.perc.250x':'BAM_PERC_250X',
                                'bam.perc.500x':'BAM_PERC_500X','bam.perc.1000x':'BAM_PERC_1000X','bam.perc.2000x':'BAM_PERC_2000X',
                                'bam.mean.insertsize':'BAM_MEAN_INSERTSIZE','bam.med.insertsize':'BAM_MED_INSERTSIZE',
                                'bam.sd.insertsize':'BAM_SD_INSERTSIZE','bam.min.insertsize':'BAM_MIN_INSERTSIZE',
                                'bam.max.insertsize':'BAM_MAX_INSERTSIZE'}

        mgi_col_renamed = {'sample':'SAMPLE','cons.per.N':'PERC_N','cons.len':'CONS_LEN','cons.perc.GC':'PERC_GC',
                                'fq.trim.pass':'FASTQ_SIZE_PASS','bam.primertrim.pass':'BAM_PRIMERTRIM_PASS','bam.perc.align':'BAM_PERC_ALIGN',
                                'bam.mean.cov':'BAM_MEAN_COV','bam.med.cov':'BAM_MED_COV','bam.max-min/mean.cov':'BAM_MAXMIN_MEANCOV',
                                'bam.perc.50x':'BAM_PERC_50X','bam.perc.100x':'BAM_PERC_100X','bam.perc.250x':'BAM_PERC_250X',
                                'bam.perc.500x':'BAM_PERC_500X','bam.perc.1000x':'BAM_PERC_1000X','bam.perc.2000x':'BAM_PERC_2000X',
                                'bam.mean.insertsize':'BAM_MEAN_INSERTSIZE','bam.med.insertsize':'BAM_MED_INSERTSIZE',
                                'bam.sd.insertsize':'BAM_SD_INSERTSIZE','bam.min.insertsize':'BAM_MIN_INSERTSIZE',
                                'bam.max.insertsize':'BAM_MAX_INSERTSIZE'}

        self.renamed_columns = {'mgi':mgi_col_renamed,'nanopore':nanopore_col_renamed,'illumina':illumina_col_renamed}
        
    def AddNewMetric(self,full_metric_path,date,tech):
        new_metric = Metric(self.sample,full_metric_path,date,tech,self.renamed_columns[tech])
        self.metrics_list.append(new_metric)

    def GetMetrics(self):
        return self.metrics_list

class Metric:
    def __init__(self,sample,path,date,tech,renamed_columns):
        self.sample = sample
        self.path = path
        self.date = date
        self.tech = tech
        self.renamed_columns_map = renamed_columns

        self.CreatePdDf()
        self.GetSamplePdDf()

    def ExtractMetrics(self):
        self.seq_len = int(self.pd_df['CONS_LEN'].values[0])

    def GetSeqLen(self):
        return(self.seq_len)

    def CreatePdDf(self):
        self.pd_df = pd.read_csv(self.path)
        self.RenameColumns()

    def RenameColumns(self):
        self.pd_df = self.pd_df.rename(columns=self.renamed_columns_map)

    def GetSamplePdDf(self):
        self.pd_df = self.pd_df.loc[self.pd_df['SAMPLE'] == self.sample.GetSampleName(),:]

    def GetSampleName(self):
        return(self.sample.GetSampleName())

    def GetTech(self):
        return(self.tech)

    def GetDate(self):
        return(self.date)

    def GetStrDate(self):
        english_day = self.date.strftime("%A")
        french_day = Utils.GetFrenchDay(english_day)

        return(self.date.strftime(french_day + " le %d %B %Y"))

    def GetPdDf(self):
        return(self.pd_df)


class Utils():
    
    @staticmethod
    def GetFrenchDay(english_day):
        day_map = {"Monday": 'Lundi', 'Tuesday': 'Mardi', 'Wednesday': 'Mercredi', 'Thursday': 'Jeudi', 'Friday': 'Vendredi', 'Saturday': 'Samedi', 'Sunday': 'Dimanche'}

        return(day_map[english_day])

    @staticmethod
    def ParseSampleDir(dir_name):
        
        pattern_obj = re.compile(r'^\S+\.(\S+)\.(\S+)\.(\d{8})$')
        try:
            search_obj = pattern_obj.search(dir_name)
            spec_name = search_obj.group(1)
            tech = search_obj.group(2)
            date = search_obj.group(3)
        except:
          print("Erreur de parsing pour ", dir_name)
          return None

        return (spec_name,tech,date)


    @staticmethod
    def GetObjectByName(obj_name,obj_list):
        for obj in obj_list:
            if obj.GetSampleName() == obj_name:
                return obj

class MarkdownWriter(object):
    author = "Sandrine Moreira"
    current_date = datetime.today().strftime('%Y-%m-%d') 
    out_type = "html_document"

    if _DEBUG:
        out_dir = "/data/Users/Eric/Covid19/TestRapportSequence/MARKDOWN_HTML_REPORT"
    else:
        out_dir = ""
 
    @classmethod 
    def GetHeader(cls,sample_name):
        key_val = {"titre" : "Rapport pour " + sample_name,"auteur" : cls.author, "date": cls.current_date, "out_type" : cls.out_type}
        header = "---\ntitle: \"{0[titre]}\"\nauthor: \"{0[auteur]}\"\ndate: \"{0[date]}\"\noutput: {0[out_type]}\n---\n".format(key_val)
        return(header)

    @classmethod
    def GetInfo(cls,metric):

        sample_name_line = "### Identifiant de l'échantillon : {0} ".format("<span style=\"color:blue\">" + metric.GetSampleName() + "<span>")
        sample_date_line = "### Date de collecte: {0} ".format("<span style=\"color:blue\">" + " a venir " + "</span>")
        seq_tech_line = "### Technologie de séquencage: {0} ".format("<span style=\"color:blue\">" + metric.GetTech() + "</span>")
        seq_date = "### Date de séquencage: {0} ".format("<span style=\"color:blue\">" + metric.GetStrDate() + "</span>")
        seq_len = "### Longueur de séquence: {0} ".format("<span style=\"color:blue\">" + str(metric.GetSeqLen()) + "</span>")
        clade_line = "### Clade: {0} ".format("<span style=\"color:blue\">" + "a venir" + "</span>")

        info = "{name}\n{samp_date}\n{tech}\n{seq_date}\n{seq_len}\n{clade}".format(name=sample_name_line,samp_date=sample_date_line,tech=seq_tech_line,seq_date=seq_date,seq_len=seq_len,clade=clade_line)

        return(info)

    @classmethod
    def GetVariantTable(cls,metric):

        title_line = "### Variants\n"
        header = "| Position | Variation nucléotidique |Variation en acide aminé | Qualité |"
        separator = "| :----: | :----: | :----: | :----: |"
        test_line = "| position | nuc_refPOSITIONnuc_var |aa_refPOSITIONaa_var | valeur de qualité |"

        table = "{title}\n{header}\n{sep}\n{test}".format(title=title_line,header=header,sep=separator,test=test_line)

        return(table)


    @classmethod
    def BuildSeqReport(cls,metric):
        sample_name = metric.GetSampleName()
        tech = metric.GetTech()
        seq_date = metric.GetDate()
        seq_date_str = seq_date.strftime("%Y%m%d")

        today = datetime.today()

        metric.ExtractMetrics()

        header = cls.GetHeader(sample_name)
        info = cls.GetInfo(metric)
        variant_table = cls.GetVariantTable(metric)

        rmd_file_name = sample_name + "_" + tech + "_" + seq_date_str + ".Rmd"
        rmd_out = os.path.join(cls.out_dir,rmd_file_name)
        cls.WriteRmd(rmd_out,header,info,variant_table)

        html_file_name = sample_name + "_" + tech + "_" + seq_date_str + ".html"
        html_out = os.path.join(cls.out_dir,html_file_name)

        cls.WriteHtml(html_out,rmd_out)

    @classmethod
    def WriteRmd(cls,out,header,info,variant_table):
        out_handler = open(out,'w')
        out_handler.write(header)
        out_handler.write("\n")
        out_handler.write(info)
        out_handler.write("\n")
        out_handler.write(variant_table)
        out_handler.write("\n")

        out_handler.close()

    @classmethod
    def WriteHtml(cls,html_out,rmd_in):
        os.system("R -e \"rmarkdown::render('{0}',output_file='{1}') \" > /dev/null".format(rmd_in,html_out))
        os.chmod(html_out,stat.S_IRWXU|stat.S_IRWXU|stat.S_IRWXU)
        

def BuildSeqReports(plate_manager):
    plate_manager.SetSampleList()
    sample_list = plate_manager.GetSamples()

    for samples in sample_list:
        samples.BuildSeqReport()


def Main():
    global basedir
    basedir = PlateDirManager.GetBaseDir(_DEBUG)
   
    plate_manager = PlateManager()
 
    for plate in os.listdir(basedir):
        plate_manager.AddPlate(plate)

    BuildSeqReports(plate_manager)

if __name__ == '__main__':
    Main()
