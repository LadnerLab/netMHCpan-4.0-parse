# netMHCpan-4.0-parse
Parses output produced by netMHCpan-4.0, and stores it in an easy to read tab-delimited file.

### Dependencies
Python >= 3.4

### Input
1. xls file created by netMHCpan-4.0
2. floating-point threshold determining what makes a strong  match
3. floating-point threshold determining what makes a weak  match
4. file to write output to


### Output
1. tab-delimited text file summary statistics created by netMHCpan-4.0 that is most 
easily viewed in a spreadsheet program

Example output:
```
Pep-len	Allele	Total	#SB	#WB
8	HLA-A02:01	2	0	2
9	HLA-A02:01	5	2	3
10	HLA-A02:01	5	0	5
11	HLA-A02:01	1	0	1
8	SLA-3-CDY	2	0	2
9	SLA-3-CDY	24	6	18
10	SLA-3-CDY	1	0	1
11	SLA-3-CDY	0	0	0
8	SLA-6:0101	0	0	0
9	SLA-6:0101	27	10	17
10	SLA-6:0101	3	1	2
11	SLA-6:0101	0	0	0
8	SLA-6:0102	0	0	0
9	SLA-6:0102	27	10	17
10	SLA-6:0102	3	1	2
11	SLA-6:0102	0	0	0
```

### Usage
```

usage: parse_outputs.py [-h] [-f IN_FILE] [-w WEAK] [-s STRONG] [-o OUTPUT]
                        [-v]

Parse output produced by netMHCpan-4.0

optional arguments:
  -h, --help            show this help message and exit
  -f IN_FILE, --in_file IN_FILE
                        Output file produced by netMHCpan-4.0 to parse
  -w WEAK, --weak WEAK  Floating-point threshold for a binding to be
                        considered weak. Anything less than this value will be
                        considered a weak binding if it is also greater than
                        or equal the strong threshold. (0.5, 0.2]
  -s STRONG, --strong STRONG
                        Floating-point threshold for a binding to be
                        considered strong. Anything less than this, will be
                        considered a strong binding. (0, 0.2)
  -o OUTPUT, --output OUTPUT
                        File to write tab-delimited output to
  -v, --verbose         Flag for output to be written to STDOUT as well as
                        output file
```

### Example command 
``` bash
./parse_outputs.py -f netMHC_out.xls -w 2.0 -s 0.2 -o MHC_out_parsed.txt
```
