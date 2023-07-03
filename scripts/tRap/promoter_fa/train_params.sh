#learn the parameters of the generalized extreme value (GEV) distributions from an appropriate training set.
#comstruce 2 new tRap packages including new parameters and new PWM: tRap_jaspar2022, tRap_pwmhocomoco11

#1 download data sets of all (putative) promoter sequences of lengths 200, 500, 600, 800, 1000, 2000, 5000 from UCSC Genome Browser(http://genome.ucsc.edu/)


#2 learn the parameters


#2.1 JASPAR 2022
mkdir jaspar2022;cd jaspar2022

##2.1.1 download PFM matrix from the JASPAR 2022 database
wget -c https://jaspar.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt --no-check-certificate 
grep '>' JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt > jaspar2022.name.txt
sed -i 's/>//g' jaspar2022.name.txt
awk -F '\t' '{print $1}' jaspar2022.name.txt > jaspar2022.motif.txt #841 motif

##2.1.2 split the motif list into some small files (7 promoter sequences,1 motif~0.4-0.5h)
mkdir motiffiles gev_parameters
split -l 13 jaspar2022.motif.txt -d motiffiles/motif

##2.1.3 change the code: fit-gev.R --> new-fit-gev.R
#../../tRap/R/fit-gev.R: add '#' for line3-line7
#cp ../../tRap/inst/scripts/fit-gev.R ./new-fit-gev.R
#vim new-fit-gev.R

##2.1.4 train the parameters
#run all the motifs
for i in {0..9}
do
nohup Rscript --vanilla new-fit-gev.R motif0$i &
done
for i in {10..64}
do
nohup Rscript --vanilla new-fit-gev.R motif$i &
done

##2.1.5 merge the parameters
cat gev_parameters/motif*.params.txt > summary_params.txt
grep -v 'matrix' summary_params.txt > summary.params.txt
sed -i '1i\matrix\tshape0\tshape1\tscale0\tscale1\tloc0\tloc1\tp.shape\tp.scale\tp.loc' summary.params.txt
rm summary_params.txt
cat gev_parameters/motif*.length.txt > summary_length.txt
grep -v 'matrix' summary_length.txt > summary.length.txt
sed -i '1i\matrix\tregion.size\tloc\tscale\tshape\tconvergence' summary.length.txt
rm summary_length.txt

##2.1.6 construct new tRap package including new parameters and pwm from the JASPAR 2022 database
#download tRap0.4
wget -c http://trap.molgen.mpg.de/download/TRAP_R_package/tRap_0.4.tar.gz
tar -zxf tRap_0.4.tar.gz
mv tRap tRap_jaspar2022
Rscript replace_parameters.R
cp ../../tRap/NAMESPACE tRap_jaspar2022/ # NAMESPACE file is missing
mkdir ~/miniconda3/lib/R/library/tRap_jaspar2022
R CMD INSTALL tRap_jaspar2022 --library=~/miniconda3/lib/R/library/tRap_jaspar2022


#2.2 HOCOMOCO v11
cd ../
mkdir pwm_hocomocov11;cd pwm_hocomocov11

##2.2.1 download PFM matrix(jaspar-format) from the HOCOMOCO v11 database
wget -c https://hocomoco11.autosome.org/final_bundle/hocomoco11/core/HUMAN/mono/HOCOMOCOv11_core_HUMAN_mono_jaspar_format.txt --no-check-certificate
grep '>' HOCOMOCOv11_core_HUMAN_mono_jaspar_format.txt > hocomoco11.name.txt
sed -i 's/>//g' hocomoco11.name.txt #401 motif

##2.2.2 split the motif list into some small files (7 promoter sequences,1 motif~0.4-0.5h)
mkdir motiffiles gev_parameters
split -l 14 hocomoco11.name.txt -d motiffiles/motif

##2.2.3 change the code: fit-gev.R --> new-fit-gev.R
cp ../../tRap/inst/scripts/fit-gev.R ./new-fit-gev.R
#vim new-fit-gev.R

##2.2.4 train the parameters 
#run all the motifs
for i in {0..9}
do
nohup Rscript --vanilla new-fit-gev.R motif0$i &
done
for i in {10..28}
do
nohup Rscript --vanilla new-fit-gev.R motif$i &
done

##2.2.5 merge the parameters
cat gev_parameters/motif*.params.txt > summary_params.txt
grep -v 'matrix' summary_params.txt > summary.params.txt
sed -i '1i\matrix\tshape0\tshape1\tscale0\tscale1\tloc0\tloc1\tp.shape\tp.scale\tp.loc' summary.params.txt
rm summary_params.txt
cat gev_parameters/motif*.length.txt > summary_length.txt
grep -v 'matrix' summary_length.txt > summary.length.txt
sed -i '1i\matrix\tregion.size\tloc\tscale\tshape\tconvergence' summary.length.txt
rm summary_length.txt

##2.2.6 construct new tRap package including new parameters and pwm from the JASPAR 2022 database
#download tRap0.4
wget -c http://trap.molgen.mpg.de/download/TRAP_R_package/tRap_0.4.tar.gz
tar -zxf tRap_0.4.tar.gz
mv tRap tRap_pwmhocomoco11
Rscript replace_parameters.R
cp ../../tRap/NAMESPACE tRap_pwmhocomoco11/ # NAMESPACE file is missing
mkdir ~/miniconda3/lib/R/library/tRap_pwmhocomoco11
R CMD INSTALL tRap_pwmhocomoco11 --library=~/miniconda3/lib/R/library/tRap_pwmhocomoco11
