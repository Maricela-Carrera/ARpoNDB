#  This program uses the phyper and dhyper functions from the R package 
#  version 4.2.2 to evaluate the enrichment of a particular feature within
#  a selected subset, relative to what would be expected by chance from a
#  background population. It requires four input values:
#
#     N: The total number of elements in the population
#        (e.g., all genes in Gammaproteobacteria)
#
#     H: The number of elements in the population with the feature of interest
#        interest (i.e., "successes")
#        (e.g., all genes regulated by sigma in Gammaproteobacteria)
#
#     n: The number of elements in the selected subset
#        (e.g., genes encoding proteins belonging to COG2710)
#
#     h: The number of elements in the subset with the feature of interest
#        (e.g., COG2710 genes regulated by sigma54 in Gammaproteobacteria)
#
#  The program reads a flat file named input_data.txt, which contains these 
#  values separated by tabs.
   
use Statistics::R ;
my $R = Statistics::R->new();

	my $out1 = $R->run(
	q`t <- read.csv("archivo",header = FALSE, sep = " ")`, 
	q`casos <- nrow(t)`,
        q`for (i in 1:casos)
		{ 
		d = t[i,1] - t[i,2];
 		ph = phyper(t[i,4], t[i,2], d, t[i,3]);
		dh = dhyper(t[i,4], t[i,2], d, t[i,3]);
		p = 1 - ph + dh;
		print(paste(ph,p))
		}`,
        );
print "$out1\n";
