# coding=utf-8
"""
Eric Fournier 2020-08-12

"""


import mysql.connector

class MySQLcovid19:
    host = 'localhost'
    user = 'root'
    password = 'lspq2019'
    database = 'TestCovid19v7'
    connection = None

    @classmethod
    def SetConnection(cls):
        cls.connection =  mysql.connector.connect(host=cls.host,user=cls.user,password=cls.password,database=cls.database)
  
    @classmethod
    def GetConnection(cls):
        return cls.connection

    @classmethod
    def GetCursor(cls):
        return cls.GetConnection().cursor()

    @classmethod
    def Commit(cls):
        cls.connection.commit()


class MySQLcovid19Selector:

    @classmethod
    def GetSampleDate(cls,cursor,sample_name):

        #j créé un repertoire bidon pour /data/Users/Eric/Covid19/TestRapportSequence/PLATE_DIR/LSPQ_002/LSPQ_002.L00266938.illumina.20200624 qui est present dans MySQL
        cursor.execute("select DATE_PRELEV_HOPITAL from Prelevements where GENOME_QUEBEC_REQUETE  like '%{0}%'".format(sample_name))
        #date_prelev = cursor.fetchone()
        date_prelev = cursor.fetchall()

        if len(date_prelev) == 1 :
            date_prelev = date_prelev[0][0]
            mois = Utils.GetFrenchMonth(date_prelev.strftime("%B"))
            return(date_prelev.strftime("%d {0} %Y".format(mois)))

        return("indéterminé")

class Utils():

    @classmethod
    def GetFrenchMonth(cls,english_month):
        month_map = {'January':'Janvier','February':'Février','March':'Mars','April':'Avril','May':'Mai','June':'Juin','Juillet':'July','August':'Août',
                     'September':'Septembre','October':'Octobre','November':'Novembre','December':'Décembre'} 

        return(month_map[english_month])
