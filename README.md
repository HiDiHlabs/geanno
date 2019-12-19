# pyanno: A package for region based annotation

pyanno is a python package that - given a list of genomic intervals - annotates these intervals to arbitrary genomic intervals of interest, as well as genes. The list of genomic intervals to be annotated, as well as the database interval regions of interest have to be encoded in a bed-like format.

## Example

### Base file
The base file of intervals, that will be annotated against a database of bed like files, that contain genomic intervals of interest has to be itself bed-like. It suffuces that the first three lines follow bed conventions, i.e. 
* 1st column: Chromosome
* 2nd column: Start position (0-based)
* 3rd column: End position (1-based)

The base interval entries can in addition contain an arbitrary number of additional columns.

```
# base.bed
4       128887787       128887839
4       188862197       188862251
4       185746125       185746231
```

### Database file
The datbase file is a tab separated file, containing information about the genomic regions of interest, that shall be annotated to the base file. The database file contains the following information:
* **FILENAME**: Absolute path to the file 
* **REGION.TYPE**: E.g. protein.coding.genes, Enhancers, ...
* **SOURCE**: E.g., Cell type from which regions are derived
* **ANNOTATION.BY**: SOURCE | NAME
* **MAX.DISTANCE**: Maximal distance between base and database intervall, such that database intervall is anotated to base intervall.
* **DISTANCE.TO**: The location to which the distance shall be computed. Can be START | END | MID | REGION
* **N.HITS**: Can be either of ALL | CLOSEST
* **NAME.COL**: If ANNOTATION.BY == NAME, then you can define the column (0-based) in which the name is stored. If NAME.COL == NA, then it is assumed, that the 4th column contains the name.

The first line of the tsv file must contain the above **bold** column identifiers!

```
# database.tsv
FILENAME        REGION.TYPE     SOURCE  ANNOTATION.BY   MAX.DISTANCE    DISTANCE.TO     N.HITS  NAME.COL
E045_15_coreMarks_dense_7_Enh.bed    Enhancer.Roadmap        E045    SOURCE  0       REGION  CLOSEST NA
E036_15_coreMarks_dense_7_Enh.bed    Enhancer.Roadmap        E036    SOURCE  0       REGION  CLOSEST NA
protein_coding_genes.bed gencode19.protein.coding.TSS    gencode.v19     NAME    200000  START   CLOSEST NA
enhancer_promoter_links_neutrophils.bed      PCHiC.neutrophils       neutrophils     NAME    10000   REGION  ALL     7
```

#### database bed like files

```
# E045_15_coreMarks_dense_7_Enh.bed
10      179800  180000  7_Enh   0       .       179800  180000  255,255,0
10      182800  183200  7_Enh   0       .       182800  183200  255,255,0
10      265600  266000  7_Enh   0       .       265600  266000  255,255,0
.
.
.
```

```
# E036_15_coreMarks_dense_7_Enh.bed
10      132200  132400  7_Enh   0       .       132200  132400  255,255,0
10      133400  133600  7_Enh   0       .       133400  133600  255,255,0
10      152600  153000  7_Enh   0       .       152600  153000  255,255,0
.
.
.
```

```
# protein_coding_genes.bed
#chrom  start   end     ensembl.id      score   strand  hugo.name
1       69091   70008   ENSG00000186092.4       NA      +       OR4F5
1       134901  139379  ENSG00000237683.5       NA      -       AL627309.1
1       367640  368634  ENSG00000235249.1       NA      +       OR4F29
.
.
.
```

```
# enhancer_promoter_links_neutrophils.bed
#oeChr  oeStart oeEnd   oeName  baitChr baitStart       baitEnd baitName
1       1150970 1156235 .       1       850619  874081  AL645608.1;RP11-54O7.3;SAMD11
1       1000704 1005126 .       1       903641  927394  C1orf170;PLEKHN1
1       1150970 1156235 .       1       903641  927394  C1orf170;PLEKHN1
.
.
.
```

```python
import pyanno
import pandas as pnd

database_filename = "database.tsv"
base_filename = "base.bed"
results_filename = "results.bed"

# Create a new GenomicRegionAnnotator instance
gra = pyanno.Annotator.GenomicRegionAnnotator()

# load base
gra.load_base_from_file(base_filename)

# load database
gra.load_database_from_file(database_filename)

# Annotate base against all database genomic region files
gra.annotate()

# Retrieve annotated base intervals as pandas.DataFrame instance
annotated_base_df = gra.get_base()

# Write annotated base intervals to disk
annotated_base_df.to_csv("/home/biegm/tmp/annotated_dmrs.csv", sep="\t", index=False)
```

## Help on module pyanno.Annotator in pyanno:

### class GenomicRegionAnnotator()


#### \_\_init\_\_(self)
```python
'''
	Standard Constructor. Creates an empty GenomicRegionAnnotator.
    
	args: None
    
	kwargs: None
'''
```

#### load_base_from_dataframe(self, base_dataframe)
```python
'''
    Function that loads base dataframe, that will be annotated against
    annotation database."
    
    args:
            base_dataframe: pandas.DataFrame
                    First three columns must be bed-like, i.e.
                    containing chromosome, start-, and end-
                    position. Must contain a header.
'''
```

#### load_base_from_file(self, base_filename)
```python
'''
    Function that loads base file, that will be annotated against
    annotation database.
    
    args:
            base_filename: string
                    First three columns must be bed-like, i.e.
                    containing chromosome, start-, and end-
                    position. Must contain a header.
'''
```

#### load_database_from_dataframe(self, database_dataframe)
```python
'''
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
    DISTANCE.TO: The location to which the distance shall
    be computed. Can be START | END | MID | REGION
    N.HITS: Can be either of ALL | CLOSEST
    NAME.COL: If ANNOTATION.BY == NAME, then you can define the
    column (0-based) in which the name is stored. If NAME.COL ==
    NA, then it is assumed, that the 4th column contains the name.
    
    args:
            database_dataframe: pandas.DataFrame
                    DataFrame that contains the database.
'''
```

#### load_database_from_file(self, database_filename)
```python
'''
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
    DISTANCE.TO: The location to which the distance shall
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
'''
```

#### annotate(self)
```python
'''
    Method, that annotates the base region table against the ROI tables
    in the database.
'''
```


#### get_base(self)
```python
'''
    Method that return self.__base
'''
```

#### print_base(self)
```python
'''
    Method that prints base.
'''
```

#### print_database(self)
```python
'''
    Method that prints the database.
'''
```
