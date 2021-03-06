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
- nom de mois en francais
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
from Covid19DB import MySQLcovid19,MySQLcovid19Selector


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
            vcf_file = glob.glob(full_metric_path + '/*major.vcf')[0]
            self.metric_manager.AddNewMetric(metric_file,vcf_file,date,tech) 
            
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


class VCF:
    def __init__(self,metric,vcf_path):
        self.metric = metric
        self.vcf_path = vcf_path

        self.var_dict = {}

    def ExtractNucVariants(self):
        vcf_handler = open(self.vcf_path)
        
        for line in vcf_handler:
            if not line.startswith('#'):
                var_info = line.split('\t')
                var_pos = var_info[1]
                ref_nuc = var_info[3]
                var_nuc = var_info[4]
                var_qual = var_info[6]
                self.var_dict[int(var_pos)] = {'ref_nuc':ref_nuc,'var_nuc':var_nuc,'var_qual':var_qual}

        vcf_handler.close()

    def GetVar(self):
        return(self.var_dict)


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
        
    def AddNewMetric(self,full_metric_path,full_vcf_path,date,tech):
        new_metric = Metric(self.sample,full_metric_path,full_vcf_path,date,tech,self.renamed_columns[tech])
        self.metrics_list.append(new_metric)

    def GetMetrics(self):
        return self.metrics_list

class Metric:
    def __init__(self,sample,metric_path,vcf_path,date,tech,renamed_columns):
        self.sample = sample
        self.metric_path = metric_path
        self.vcf = VCF(self,vcf_path)
        self.date = date
        self.tech = tech
        self.renamed_columns_map = renamed_columns

        self.CreatePdDf()
        self.GetSamplePdDf()


    def ExtractMetrics(self):
        self.seq_len = int(self.pd_df['CONS_LEN'].values[0])
        self.perc_N = self.pd_df['PERC_N'].values[0]
        self.perc_GC = self.pd_df['PERC_GC'].values[0]
        self.mean_cov = self.pd_df['BAM_MEAN_COV'].values[0]
        self._50x_cov = self.pd_df['BAM_PERC_50X'].values[0] 
        self._100x_cov = self.pd_df['BAM_PERC_100X'].values[0]
        self._250x_cov = self.pd_df['BAM_PERC_250X'].values[0]
        self._500x_cov = self.pd_df['BAM_PERC_500X'].values[0]
        self._1000x_cov = self.pd_df['BAM_PERC_1000X'].values[0]
        
    def ExtractNucVariants(self):
        self.vcf.ExtractNucVariants()

    def Get1000xCov(self):
        return(self._1000x_cov)

    def Get500xCov(self):
        return(self._500x_cov)

    def Get250xCov(self):
        return(self._250x_cov)

    def Get100xCov(self):
        return(self._100x_cov)

    def Get50xCov(self):
        return(self._50x_cov)

    def GetMeanCov(self):
        return(self.mean_cov)

    def GetPercGC(self):
        return(self.perc_GC)

    def GetSeqLen(self):
        return(self.seq_len)

    def GetPercN(self):
        return(self.perc_N)

    def CreatePdDf(self):
        self.pd_df = pd.read_csv(self.metric_path)
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
        
        english_month = self.date.strftime("%B")
        french_month = Utils.GetFrenchMonth(english_month)

        return(self.date.strftime(french_day + " le %d {0} %Y".format(french_month)))

    def GetPdDf(self):
        return(self.pd_df)

    def GetSampleDate(self):
        sample_date = MySQLcovid19Selector.GetSampleDate(MySQLcovid19.GetCursor(),self.sample.GetSampleName())
        return(sample_date)

    @staticmethod
    def GetPercN_Quality_Score(perc_N):
        #TODO a determiner
        if perc_N < 10:
            return "Excellent"
        elif perc_N > 10 and perc_N < 20:
            return "Bon"
        else:
            return "Faible"

    @staticmethod
    def Get50xCov_Quality_Score(_50x_cov):
        #TODO a determiner
        pass 

class Utils():
    
    @staticmethod
    def GetFrenchDay(english_day):
        day_map = {"Monday": 'Lundi', 'Tuesday': 'Mardi', 'Wednesday': 'Mercredi', 'Thursday': 'Jeudi', 'Friday': 'Vendredi', 'Saturday': 'Samedi', 'Sunday': 'Dimanche'}

        return(day_map[english_day])

    @staticmethod
    def GetFrenchMonth(english_month):
        month_map = {'January':'Janvier','February':'Février','March':'Mars','April':'Avril','May':'Mai','June':'Juin','Juillet':'July','August':'Août',
                     'September':'Septembre','October':'Octobre','November':'Novembre','December':'Décembre'}

        return(month_map[english_month])


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
        sample_date_line = "### Date de collecte: {0} ".format("<span style=\"color:blue\">" + metric.GetSampleDate() + "</span>")
        seq_tech_line = "### Technologie de séquencage: {0} ".format("<span style=\"color:blue\">" + metric.GetTech() + "</span>")
        seq_date = "### Date de séquencage: {0} ".format("<span style=\"color:blue\">" + metric.GetStrDate() + "</span>")
        seq_len = "### Longueur de séquence: {0} ".format("<span style=\"color:blue\">" + str(metric.GetSeqLen()) + "</span>")
        clade_line = "### Clade: {0} ".format("<span style=\"color:blue\">" + "a venir" + "</span>")

        info = "{name}\n{samp_date}\n{tech}\n{seq_date}\n{seq_len}\n{clade}".format(name=sample_name_line,samp_date=sample_date_line,tech=seq_tech_line,seq_date=seq_date,seq_len=seq_len,clade=clade_line)

        return(info)

    @classmethod
    def GetVariantTable(cls,metric):
        def AddVarLine(var):
            return("| {0} | {1} | na | {2} |\n".format(str(var[0]),var[1]['ref_nuc']+str(var[0])+var[1]['var_nuc'],var[1]['var_qual']))

        var_list = metric.vcf.GetVar()

        title_line = "### Variants\n"
        header = "| Position | Variation nucléotidique |Variation en acide aminé | Qualité |"
        separator = "| :----: | :----: | :----: | :----: |"
        var_lines = ""

        for var in sorted(var_list.items()):
            var_lines += AddVarLine(var)

        table = "{title}\n{header}\n{sep}\n{var_lines}".format(title=title_line,header=header,sep=separator,var_lines=var_lines)

        return(table)


    @classmethod
    def GetQualityTable(cls,metric):
        
        title_line = "### Qualité\n"
        header = "| Critères | Valeur | Qualité | Seuil |"
        separator = "| :----: | :----: | :----: | :----: |"

        perc_N_line = "| Pourcentage de bases indéterminées | {perc_n} | {qual_score} | {seuil} |".format(perc_n = "%.1f" % metric.GetPercN(),qual_score = Metric.GetPercN_Quality_Score(metric.GetPercN()),seuil = "]1% - 5%[")

        perc_cov_50x =  "| Pourcentage de couverture à minimum 50X | {_50x_cov} | {qual_score} | {seuil} |".format(_50x_cov = "%.1f" % metric.Get50xCov(),qual_score = "a déterminer",seuil = "a déterminer")

        perc_cov_100x =  "| Pourcentage de couverture à minimum 100X | {_100x_cov} | {qual_score} | {seuil} |".format(_100x_cov = "%.1f" % metric.Get100xCov(),qual_score = "a déterminer",seuil = "a déterminer")

        perc_cov_250x =  "| Pourcentage de couverture à minimum 250X | {_250x_cov} | {qual_score} | {seuil} |".format(_250x_cov = "%.1f" % metric.Get250xCov(),qual_score = "a déterminer",seuil = "a déterminer")

        perc_cov_500x =  "| Pourcentage de couverture à minimum 500X | {_500x_cov} | {qual_score} | {seuil} |".format(_500x_cov = "%.1f" % metric.Get500xCov(),qual_score = "a déterminer",seuil = "a déterminer")

        perc_cov_1000x =  "| Pourcentage de couverture à minimum 1000X | {_1000x_cov} | {qual_score} | {seuil} |".format(_1000x_cov = "%.1f" % metric.Get1000xCov(),qual_score = "a déterminer",seuil = "a déterminer")

        table = "{title}\n{header}\n{sep}\n{perc_n}\n{perc_cov_50x}\n{perc_cov_100x}\n{perc_cov_250x}\n{perc_cov_500x}\n{perc_cov_1000x}".format(title=title_line,header=header,sep=separator,perc_n=perc_N_line,perc_cov_50x=perc_cov_50x,perc_cov_100x=perc_cov_100x,perc_cov_250x=perc_cov_250x,perc_cov_500x=perc_cov_500x,perc_cov_1000x=perc_cov_1000x)

        return(table) 


    @classmethod
    def BuildSeqReport(cls,metric):
        sample_name = metric.GetSampleName()
        tech = metric.GetTech()
        seq_date = metric.GetDate()
        seq_date_str = seq_date.strftime("%Y%m%d")

        today = datetime.today()

        metric.ExtractMetrics()
        metric.ExtractNucVariants()

        header = cls.GetHeader(sample_name)
        info = cls.GetInfo(metric)
        variant_table = cls.GetVariantTable(metric)
        quality_table = cls.GetQualityTable(metric)

        rmd_file_name = sample_name + "_" + tech + "_" + seq_date_str + ".Rmd"
        rmd_out = os.path.join(cls.out_dir,rmd_file_name)
        cls.WriteRmd(rmd_out,header,info,variant_table,quality_table)

        html_file_name = sample_name + "_" + tech + "_" + seq_date_str + ".html"
        html_out = os.path.join(cls.out_dir,html_file_name)

        cls.WriteHtml(html_out,rmd_out)

    @classmethod
    def WriteRmd(cls,out,header,info,variant_table,quality_table):
        out_handler = open(out,'w')
        out_handler.write(header)
        out_handler.write("\n")
        out_handler.write(info)
        out_handler.write("\n")
        out_handler.write("\n")
        out_handler.write(variant_table)
        out_handler.write("\n")
        out_handler.write("\n")
        out_handler.write(quality_table)
        out_handler.write("\n")
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

    MySQLcovid19.SetConnection()

    global basedir
    basedir = PlateDirManager.GetBaseDir(_DEBUG)
   
    plate_manager = PlateManager()
 
    for plate in os.listdir(basedir):
        plate_manager.AddPlate(plate)

    BuildSeqReports(plate_manager)

if __name__ == '__main__':
    Main()
