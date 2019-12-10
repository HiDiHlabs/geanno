import pybedtools
import pandas as pnd
from copy import deepcopy
import os

class GenomicRegionAnnotator():
	#############################
	# Constructors/ Destructors #
	#############################

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
		# SOURCE: E.g., Cell type from which regions are derived
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


	##################
	# Public methods #
	##################

	######################
	# Data loading methods

	def load_database_from_file(self, database_filename):
		'''
			Method for loading a database from a tab separated file.
			The database contains all files against which the annotation
			shall be performed. Required columns are
			FILENAME: Absolute path to the file (must be a bed like file)
			REGION.TYPE: E.g. protein.coding.genes, Enhancers, ...
			SOURCE: E.g., Cell type from which regions are derived
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
		required_fields = ["FILENAME", "REGION.TYPE", "SOURCE", "ANNOTATION.TYPE", 
					"MAX.DISTANCE", "DISTANCE.TO"]
		for required_field in required_fields:
			if(not required_field in columns):
				self.__database = None
				raise(TypeError(required_field+(" is a required column in the database "
                                                        "DataFrame but not found. Please update "
                                                        "your database DataFrame!")))

		# Check if files in __database exist
		self.__check_database()

	def load_database_from_dataframe(self, database_dataframe):
		'''
			Method for loading a database from a pandas.DataFrame.
			The database contains all files against which the annotation
			shall be performed. Required columns are
			FILENAME: Absolute path to the file 
			REGION.TYPE: E.g. protein.coding.genes, Enhancers, ...
			SOURCE: E.g., Cell type from which regions are derived
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
		required_fields = ["FILENAME", "REGION.TYPE", "SOURCE", "ANNOTATION.TYPE", 
					"MAX.DISTANCE", "DISTANCE.TO"]
		for required_field in required_fields:
			if(not required_field in columns):
				self.__database = None
				raise(TypeError(required_field+(" is a required column in the database "
                                                        "DataFrame but not found. Please update "
                                                        "your database DataFrame!")))
		
		# Check if files in __database exist
		self.__check_database()

	def load_base_from_file(self, base_filename):
		'''
			Function that loads base file, that will be annotated against
			annotation database.

			args:
				base_filename: string
					First three columns must be bed-like, i.e.
					containing chromosome, start-, and end-
					position. Must contain a header.
		'''
		self.__base = pnd.read_csv(base_filename, sep="\t")
		# Check if __base contains header and is bed-like
		self.__check_base()

		self.__base.index = [ "_".join([ str(e) for e in r.iloc[:3] ]) for i, r in self.__base.iterrows() ]

		# Create pybedtools.BedTool pbject from self.__base
		self.__base_bed = self.__create_bed4(self.__base)

	def load_base_from_dataframe(self, base_dataframe):
		'''
			Function that loads base dataframe, that will be annotated against
			annotation database.

			args:
				base_dataframe: pandas.DataFrame
					First three columns must be bed-like, i.e.
					containing chromosome, start-, and end-
					position. Must contain a header.
		'''
		self.__base = deepcopy(base_dataframe)
		# Check if __base contains header and is bed-like
		self.__check_base()

		self.__base.index = [ "_".join([ str(e) for e in r.iloc[:3] ]) for i, r in self.__base.iterrows() ]
		
		# Create pybedtools.BedTool pbject from self.__base
		self.__base_bed = self.__create_bed4(self.__base)

	###############
	# Print methods

	def print_database(self):
		'''
			Method that prints the database.
		'''
		print(self.__database)

	def print_base(self):
		'''
			Method that prints base.
		'''
		print(self.__base)


	###################
	# Private Methods #
	###################

	def __check_database(self):
		'''
			Checks if files in database exist!
		'''
		if(self.__database is None):
			print("No annotation database defined yet!")
		else:
			for index, row in self.__database.iterrows():
				filename = row["FILENAME"]
				if(not(os.path.exists(filename))):
					raise(IOError(filename+" does not exist!"))

	def __check_base(self):
		'''
			Checks if base, that will be annotated is bed-like, and
			contains a header
		'''
		if(not self.__base.columns[0][0] == "#"):
			raise(ValueError(("Base table does not contain a "
					"valid haeder! Header has to start with \"#\"")))
		for index, row in self.__base.iterrows():
			if(not(type(row.iloc[1]) == int and type(row.iloc[2]) == int)):
				raise(TypeError(("Base table does not seem to be bed-like. "
						"Second and third columns must be integers.")))
			elif(not(row.iloc[2] > row.iloc[1])):
				raise(TypeError(("Base table does not seem to be bed-like. "
						"Second column must be smaller or equal to third "
						"column.")))

	def __create_bed4(self, df):
		'''
			Method that creates a bed4 pybedtools.BedTool object
			from df. Columns are: 1. Chromosome, 2. Start,
			3. End, 4. Name (<chrom>_<start>_<end>)

			args:
				df: pnd.DataFrame
		'''
		bed_list = []
		
		for index, row in df.iterrows():
			bed_list += [ "\t".join([ str(e) for e in row.iloc[:3] ]+
					[ "_".join([ str(e) for e in row.iloc[:3] ]) ]) ]

		return pybedtools.BedTool("\n".join(bed_list), from_string=True)

	def __anno_done(self, region_type, source, annotation_type):
		'''
			Method that checks if annotation is already done for
			region_type, source, annotation_type combo.
		'''
		if(not region_type in set(self.__base.index)):
			return False
		elif( annotation_type in set(self.__base.index) ):
			sources = set(sum([], [ e.split(";") for e in self.__base.loc[:, region_type] ]))
			if(not source in sources):
				return False
		else:
			return True	
