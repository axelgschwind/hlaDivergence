# hlaDivergence
This tool calculates protein distances between pairs of human class I and class II HLA alleles.

## Usage
HLA pairs can be passed via commandline in four-digit nomenclature. e.g

    > python main.py HLA-A*01:01 HLA-A*02:01

HLA pairs can also be processed from a batch file which contains the HLA pairs tab-separated per row:

    > python main.py --batch data/grouped_common_alleles.tsv
    

    
## Output

hlaDivergence prints the result to the standard output. Every line reflects the result of an HLA pair, e.g.
| #ALLELE1 | ALLELE2 | GRANTHAM | PDISTANCE | SANDBERG | DAYHOFF | JTT |
| -------- | ------- | -------- | --------- | -------- | ------- | --- |
| HLA-A*01:01 | HLA-A*02:01 | 10.50276 | 0.1326 | 0.71498 | 0.14137 | 0.12 |

## Metrics
hlaDivergence calculates the protein distance with five different metrics:
* p-Distance
* Grantham distance
* Sandberg distance
* Dayhoff distance
* JTT distance
