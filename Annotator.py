import pybedtools
import pandas as pnd
from copy import deepcopy

class GenomicRegionAnnotator():
	def __init__(self):
	'''
		Standard Constructor. Creates an empty GenomicRegionAnnotator.

		args: None

		kwargs: None
	'''
		# Set database against which to None. Will be pandas.DataFrame
		# object containing the following columns: 
		# FILENAME: Absolute path to the file 
		# REGION.TYPE: E.g. protein.coding.genes, Enhancers, ...
		# FILE.ID: E.g., Cell type from which regions are derived
		# ANNOTATION.TYPE: overlap | distance
		# MAX.DISTANCE: If ANNOTATION.TYPE is distance, then the max
		# distance has to be defined
		# DISTANCE.TO: If ANNOTATION.TYPE is distance, then it has to
		# be defined what the location is to which the distance shall
		# be computed. Can be START | END | MID
		self.__database = None

		# Set regions that shall be annotated to None.
		# List will be filled with filename, and pybedtools.BedTool
		# object
		self.__base = None
		
	def load_database_from_file(self, database_filename):
		'''
			Method for loading a database from a tab separated file.
			The database contains all files against which the annotation
			shall be performed. Required columns are
			FILENAME: Absolute path to the file 
			REGION.TYPE: E.g. protein.coding.genes, Enhancers, ...
			FILE.ID: E.g., Cell type from which regions are derived
			ANNOTATION.TYPE: overlap | distance
			MAX.DISTANCE: If ANNOTATION.TYPE is distance, then the max
			distance has to be defined
			DISTANCE.TO: If ANNOTATION.TYPE is distance, then it has to
			be defined what the location is to which the distance shall
			be computed. Can be START | END | MID

			args:
				database_filename: string
					Absolute pathe to the tab separated file
					containing information about the files that
					shall be annotated to the base file.
		'''
		self.__database = pnd.read_csv(database_filename, sep="\t")

		# Check if necessary fields are defined
		columns = set(self.__database.columns)
		required_fields = ["FILENAME", "REGION.TYPE", "FILE.ID", "ANNOTATION.TYPE", "MAX.DISTANCE", "DISTANCE.TO"]
		for required_field in required_fields:
			if(not required_field in columns):
				self.__database = None
				print(required_field+" is a required column in the database file but not found in "+datbase_filename+". Please update your database file!")
				break

	def load_database_from_dataframe(self, database_dataframe):
		'''
			Method for loading a database from a pandas.DataFrame.
			The database contains all files against which the annotation
			shall be performed. Required columns are
			FILENAME: Absolute path to the file 
			REGION.TYPE: E.g. protein.coding.genes, Enhancers, ...
			FILE.ID: E.g., Cell type from which regions are derived
			ANNOTATION.TYPE: overlap | distance
			MAX.DISTANCE: If ANNOTATION.TYPE is distance, then the max
			distance has to be defined
			DISTANCE.TO: If ANNOTATION.TYPE is distance, then it has to
			be defined what the location is to which the distance shall
			be computed. Can be START | END | MID

			args:
				database_dataframe: pandas.DataFrame
					DataFrame that contains the database.
		'''
		self.__database = deepcopy(database_dataframe)

		# Check if necessary fields are defined
		columns = set(self.__database.columns)
		required_fields = ["FILENAME", "REGION.TYPE", "FILE.ID", "ANNOTATION.TYPE", "MAX.DISTANCE", "DISTANCE.TO"]
		for required_field in required_fields:
			if(not required_field in columns):
				self.__database = None
				print(required_field+" is a required column in the database DataFrame but not found. Please update your database DataFrame!")
				break

	def print_database(self):
		print(self.__database)
