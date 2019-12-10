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

		# Set pybedtools.BedTool object, that contains data of 
		# self.__base as bed4 BedTool object to None
		self.__base_bed = None


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

	####################
	# Annotation methods
	def annotate(self):
		'''
			Method, that annotats the base region table against the ROI tables
			in the database.
		'''
		# Check if all necessary objects are defined
		if(self.__base is None):
			raise(NameError(("Base regions are not defined! Please define them "
					"using either of load_base_from_file or "
					"load_base_from_dataframe method.")))
		elif(self.__database is None):
			raise(NameError(("Database regions are not defined! Please define them "
					"using either of load_database_from_file or "
					"load_database_from_dataframe method.")))

		# Perform Annotation
		for index, row in self.__database.iterrows():
			# Check if annotation was already performed and continue if True
			if(self.__anno_done(row["REGION.TYPE"], row["SOURCE"], row["ANNOTATION.TYPE"])):
				print(("Annotation DONE for "+row["REGION.TYPE"]+" "+row["SOURCE"]+
					" "+row["ANNOTATION.TYPE"]))
				continue
			else:
				# Perform overlap based annotation
				if(row["ANNOTATION.TYPE"] == "overlap"):
					if(not row["REGION.TYPE"] in self.__base.columns):
						# Initiate new column if self.__base with "NA"
						self.__base[row["REGION.TYPE"]] = ["NA"]*len(self.__base.index)

					anno_bed = pybedtools.BedTool(row["FILENAME"])
					intersect_bed = self.__base_bed.intersect(anno_bed, wa=True, u=True)
					for e in intersect_bed:
						anno_string = self.__base.loc[e[3], row["REGION.TYPE"]]
						if(not anno_string == "NA"):
							anno_string += ";"+row["SOURCE"]
						else:
							anno_string = row["SOURCE"]
						self.__base.loc[e[3], row["REGION.TYPE"]] = anno_string

				# Perform distance based annotation
				elif(row["ANNOTATION.TYPE"] == "distance"):
					# Initiate new column if self.__base with "NA"
					self.__base[row["REGION.TYPE"]] = ["NA"]*len(self.__base.index)
					
					anno_bed = self.__create_bed6_single_base(row["FILENAME"], row["DISTANCE.TO"]).sort()
					closest_bed = self.__base_bed.sort().closest(anno_bed, d=True)

					for e in closest_bed:
						name_roi = str(e[7])
						chrom_roi = str(e[4])
						start_roi = int(e[5])
						end_roi = int(e[6])
						strand_roi = str(e[9])
						distance = int(e[-1])
						if(not chrom_roi == "none"):
							if(not distance == 0):
								if(strand_roi == "+"):
									if(start_roi > int(e[2])):
										distance = -1*distance
								elif(strand_roi == "-"):
									if(end_roi < int(e[1])):
										distance = -1*distance
						anno_string = self.__base.loc[e[3], row["REGION.TYPE"]]
						if(not(anno_string == "NA")):
							anno_string += ";"+name_roi+"("+str(distance)+")"
						else:
							anno_string = name_roi+"("+str(distance)+")"
						self.__base.loc[e[3], row["REGION.TYPE"]] = anno_string

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

	################
	# Getter methods

	def get_base(self):
		'''
			Method that return self.__base
		'''
		return deepcopy(self.__base)

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

	def __create_bed6_single_base(self, bed_filename, pos):
		'''
			Create a bed6 pybedtools.BedTool object using single base
			at start-, end- or midpoint of region.

			args:
				bed_filename: string
					Path to bed file used for creating BedTool
					object.
				pos: string
					base position relative to intervall used 
					for creating the BedTool object. Can be
					either of START | END | MID. START, and
					END are relative to strand.
		'''
		bed_list = []
		bed_file = open(bed_filename, "r")
		for line in bed_file:
			if(line[0] == "#"):
				continue
			split_line = line.rstrip().split("\t")
			chrom = split_line[0]
			start = int(split_line[1])
			end = int(split_line[2])
			name = split_line[3]
			strand = split_line[5]

			if(pos == "START"):
				if(strand == "+"):
					end = start+1
				else:
					start = end-1
			elif(pos == "END"):
				if(strand == "-"):
					start = end-1
				else:
					end = start+1
			elif(pos == "MID"):
				start = start+int((end-start)/2)
				end = start + 1
			bed_list += [ "\t".join([chrom, str(start), str(end), name, "NA", strand]) ]

		return pybedtools.BedTool("\n".join(bed_list), from_string=True)

	def __anno_done(self, region_type, source, annotation_type):
		'''
			Method that checks if annotation is already done for
			region_type, source, annotation_type combo.
		'''
		if(not region_type in set(self.__base.columns)):
			return False
		elif( annotation_type == "overlap" ):
			sources = set(sum([ e.split(";") for e in self.__base.loc[:, region_type] ], []))
			if(not source in sources):
				return False
		else:
			return True
