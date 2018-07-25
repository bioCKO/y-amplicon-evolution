Partial CNV Detection contains:

amplicon_binary_segmentation.py: Takes bedgraph coverage file and returns the maximally significant changepoint in depth of each amplicon using a modified version of the binary segmentation algorithm. Also plots amplicon depth by 100-bp windows and shows location of maximally significant changepoint.

Coverage files should be GC-corrected first if possible (see GC Correction directory). One or more GC-corrected output bedgraph files are the input of this script. Other arguments are:

-o, output path and filename for the file containing the maximally significant changepoint in depth of each amplicon and control region (default: ./amplicon_partials) and filename prefix for plots if -p is used

-s, calculate the maximally significant changepoint in depth for each control region (including the 1-Mb normalization region) and amplicon. At least one of -s or -p must be used, or the script will do nothing. When -s is used, a file will be created (default name: ./amplicon_partials) in the following format, with fields separated by tabs:

For control regions:
1. Region name, 2. Mann-Whitney U p-value of maximally significant changepoint, 3. location of maximally significant changepoint (counted in 100-bp windows), 4. mean depth to left of changepoint, 5. mean depth to right of changepoint

For amplicons:
1. Amplicon name, 2. amplicon copy name, 3. Mann-Whitney U p-value of maximally significant changepoint, 4. location of maximally significant changepoint (counted in 100-bp windows), 5. mean depth to left of changepoint, 6. mean depth to right of changepoint, 7. copy number call to left of changepoint, 8. copy number call to right of changepoint, 9+. fields 2-8 repeat for each copy of the amplicon

-p, plots the depth of each copy of amplicons in non-overlapping 100-bp windows. At least one of -s or -p must be used, or the script will do nothing. To use this option, specify which amplicon(s) to plot after the "-p" flag. Valid amplicons are IR1, IR2, IR3, IR5, P8, P7, P6, P5, P4, Blue, Teal, Green, Red, Gray, and Yellow. If -s is also used, plots will also contain a line showing the maximally significant changepoint, lines showing mean depth to the left and right of the changepoint, copy number calls to the left and right of the changepoint, and Mann-Whitney U p-value of maximally significant changepoint for each amplicon copy.

-g, a file of Y chromosome locations that are not masked by the repeat masking pipeline in the Amplicon Annotation and Repeat Masking directory. By default, this script uses the "good_chrY_amplicon_bases_w100_c11.pickle" file in that directory. (If the directory structure of the master directory isn't present, this script will not be able to locate those files, and you will have to input the location of the file using the -g flag.)

-r, a file of Y chromosome amplicon annotations (provided in the Amplicon Annotation and Repeat Masking directory as "Y_repeats.txt"). This script uses that file by default. As with -g, the directory structure must be the same as the default or else the location of this file must be specified directly using the -r flag.


usage: amplicon_binary_segmentation.py [-h] [-o OUTFILE] [-s]
                                       [-p [PLOT_DEPTHS [PLOT_DEPTHS ...]]]
                                       [-g GOOD_LOCATIONS] [-r REPEATS_FILE]
                                       [input_files [input_files ...]]

positional arguments:
  input_files           bedGraph file(s) of Y chromosome depth

optional arguments:
  -h, --help            show this help message and exit
  -o OUTFILE, --outfile OUTFILE
                        Output filename
  -s, --segmentation    Use binary segmentation to find maximum likelihood
                        breakpoint of amplicons (at least one of -s and -p
                        must be chosen)
  -p [PLOT_DEPTHS [PLOT_DEPTHS ...]], --plot_depths [PLOT_DEPTHS [PLOT_DEPTHS ...]]
                        create plots of amplicon depth (input is e.g. Blue,
                        Teal, etc.)
  -g GOOD_LOCATIONS, --good_locations GOOD_LOCATIONS
                        File of unmasked Y chromosome locations (provided in
                        Amplicon Annotation and Repeat Masking
  -r REPEATS_FILE, --repeats_file REPEATS_FILE
                        File of Y chromosome amplicon locations (provided in
                        Amplicon Annotation and Repeat Masking