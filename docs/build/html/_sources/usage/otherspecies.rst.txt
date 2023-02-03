Working with other species
__________________________

*pyTEnrich* has been developed to be used with hg19 genome. If you wish to use it with another genome, you will have to prepare the input TE database with a proper format and make new associated summaries with provided perl script.

Exemple for species danRer10 :

1. Download TE database from RepeatMasker

2. Convert \.fa\.out file with utils/convert_repeatMasker.sh
``utils/convert_repeatMasker.sh danRer10.fa.out > db/danRer10.bed``

3. Remove unwanted repeats : 
``grep -v 'Low_complexity\|Satellite\|Simple_repeat\|Unknown' db/danRer10.bed > db/danRer10_TEs.bed``

4. Use utils/make_ref_TE.pl to make files needed by pyTEnrich
    INDEX is the col number of the desired feature (subfam/fam)

*class/families should always be in field 6 (0-based)*
``utils/make_ref_TE.pl --file db/danRer10_TEs.bed --index 6 > db/danRer10_Fam_ref_TE.txt``

*subfam should always be in field 7 (0-based)*
``utils/make_ref_TE.pl --file db/danRer10_TEs.bed --index 7 > db/danRer10_Subfam_ref_TE.txt``


