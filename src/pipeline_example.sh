#build graph on initial backbone, run DBS
python AccuVIR_main.py -r ../test_5_7/nano_50x_4k/4k_50x_ec.fa -b ../test_5_7/HXB2.fasta -m 1 --beamwidth 500
#build graph on DBS longest path, run DBS and sampling
mv ../test_5_7/nano_50x_4k/4k_50x_ec.fa_ON_HXB2_DBS_500_longest.fa ../test_5_7/nano_50x_4k/longest.fa
python AccuVIR_main.py -r ../test_5_7/nano_50x_4k/4k_50x_ec.fa -b ../test_5_7/nano_50x_4k/longest.fa -m 3
#run GM on the merged paths
#run MRR to output final path
python AccuVIR_MRR.py -r ../test_5_7/nano_50x_4k/4k_50x_ec.fa_ON_longest_merge.fa
