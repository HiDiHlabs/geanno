# pyanno: A package for region based annotation

pyanno is a python package that - given a list of genomic intervals - annotates these intervals to arbitrary genomic intervals of interest, as well as genes. The list of genomic intervals to be annotated, as well as the database interval regions of interest have to be encoded in a bed-like format.

## Help on module pyanno.Annotator in pyanno:

### class GenomicRegionAnnotator()


#### \_\_init\_\_(self)

_Standard Constructor. Creates an empty GenomicRegionAnnotator._
    
	+ args: None
    
	+ kwargs: None

```python
load_base_from_dataframe(self, base_dataframe)
    "Function that loads base dataframe, that will be annotated against
    annotation database."
    
    args:
            base_dataframe: pandas.DataFrame
                    "First three columns must be bed-like, i.e.
                    containing chromosome, start-, and end-
                    position. Must contain a header."
```

```python
load_base_from_file(self, base_filename)
    Function that loads base file, that will be annotated against
    annotation database.
    
    args:
            base_filename: string
                    First three columns must be bed-like, i.e.
                    containing chromosome, start-, and end-
                    position. Must contain a header.
```

```python
load_database_from_dataframe(self, database_dataframe)
    Method for loading a database from a pandas.DataFrame.
    The database contains all files against which the annotation
    shall be performed. Required columns are
    FILENAME: Absolute path to the file 
    REGION.TYPE: E.g. protein.coding.genes, Enhancers, ...
    SOURCE: E.g., Cell type from which regions are derived
    ANNOTATION.BY: SOURCE | NAME
    MAX.DISTANCE: Maximal distance between base and database
    intervall, such that database intervall is anotated to base
    intervall.
    DISTANCE.TO: If ANNOTATION.TYPE is distance, then it has to
    be defined what the location is to which the distance shall
    be computed. Can be START | END | MID | REGION
    N.HITS: Can be either of ALL | CLOSEST
    NAME.COL: If ANNOTATION.BY == NAME, then you can define the
    column (0-based) in which the name is stored. If NAME.COL ==
    NA, then it is assumed, that the 4th column contains the name.
    
    args:
            database_dataframe: pandas.DataFrame
                    DataFrame that contains the database.
```

```python
load_database_from_file(self, database_filename)
    Method for loading a database from a tab separated file.
    The database contains all files against which the annotation
    shall be performed. Required columns are
    FILENAME: Absolute path to the file (must be a bed like file)
    REGION.TYPE: E.g. protein.coding.genes, Enhancers, ...
    SOURCE: E.g., Cell type from which regions are derived
    ANNOTATION.BY: SOURCE | NAME
    MAX.DISTANCE: Maximal distance between base and database
    intervall, such that database intervall is anotated to base
    intervall.
    DISTANCE.TO: If ANNOTATION.TYPE is distance, then it has to
    be defined what the location is to which the distance shall
    be computed. Can be START | END | MID | REGION
    N.HITS: Can be either of ALL | CLOSEST
    NAME.COL: If ANNOTATION.BY == NAME, then you can define the
    column (0-based) in which the name is stored. If NAME.COL ==
    NA, then it is assumed, that the 4th column contains the name.
    
    args:
            database_filename: string
                    Absolute pathe to the tab separated file
                    containing information about the files that
                    shall be annotated to the base file.
```

```python
annotate(self)
    Method, that annotates the base region table against the ROI tables
    in the database.
```


```python
get_base(self)
    Method that return self.__base
```

```python
print_base(self)
    Method that prints base.
```

```python
print_database(self)
    Method that prints the database.
```
