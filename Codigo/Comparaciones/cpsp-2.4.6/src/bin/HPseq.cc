
#include <iostream>
#include <vector>
#include <map>
#include <string>

#include <biu/OptionParser.hh>

#include "version.hh"

/* 
 * HPseq converts an amino acid sequence into an HP sequence using a 
 * hydrophobicity plot.
 * 
 */

////////////////////////////////////////////////////////////////////////////

//! Inits the allowed parameters for HPconvert.
void initAllowedArguments(	biu::OptionMap & allowedArgs,
							std::string &infoText ) ;
									
////////////////////////////////////////////////////////////////////////////

 //! number of amino acid one letter codes in the hydrophibicity tables
const size_t aminoAcidNumber = 20;

 //! the position coding used in the hydrophicity scale tables below
const char aminoAcidCode[aminoAcidNumber] = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'};

const double KyteDoolittle[aminoAcidNumber] = { 1.80, 2.50, -3.50, -3.50, 2.80, -0.40, -3.20, 4.50, -3.90, 3.80, 1.90, -3.50, -1.60, -3.50, -4.50, -0.80, -0.70, 4.20, -0.90, -1.30 };

const double HoppWoodsNegated[aminoAcidNumber] = { 0.50, 1.00, -3.00, -3.00, 2.50, 0.00, 0.50, 1.80, -3.00, 1.80, 1.30, -0.20, 0.00, -0.20, -3.00, -0.30, 0.40, 1.50, 3.40, 2.30 };

const double Cornette[aminoAcidNumber] = { 0.20, 4.10, -3.10, -1.80, 4.40, 0.00, 0.50, 4.80, -3.10, 5.70, 4.20, -0.50, -2.20, -2.80, 1.40, -0.50, -1.90, 4.70, 1.00, 3.20 };

const double Eisenberg[aminoAcidNumber] = { 0.62, 0.29, -0.90, -0.74, 1.19, 0.48, -0.40, 1.38, -1.50, 1.06, 0.64, -0.78, 0.12, -0.85, -2.53, -0.18, -0.05, 1.08, 0.81, 0.26 };

const double Rose[aminoAcidNumber] = { 0.74, 0.91, 0.62, 0.62, 0.88, 0.72, 0.78, 0.88, 0.52, 0.85, 0.85, 0.63, 0.64, 0.62, 0.64, 0.66, 0.70, 0.86, 0.85, 0.76 };

const double Janin[aminoAcidNumber] = { 0.30, 0.90, -0.60, -0.70, 0.50, 0.30, -0.10, 0.70, -1.80, 0.50, 0.40, -0.50, -0.30, -0.70, -1.40, -0.10, -0.20, 0.60, 0.30, -0.40 };

const double EngelmanGES[aminoAcidNumber] = { 1.60, 2.00, -9.20, -8.20, 3.70, 1.00, -3.00, 3.10, -8.80, 2.80, 3.40, -4.80, -0.20, -4.10, -12.3, 0.60, 1.20, 2.60, 1.90, -0.70 };
									
////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv )
{
	  // init allowed parameters
	biu::OptionMap allowedArgs;
	std::string infoText;
	initAllowedArguments( allowedArgs, infoText );
	
	  // parse calling parameters
	biu::COptionParser opts( allowedArgs, argc, argv, infoText );
	if (!opts.noErrors()) {
		  // arguments not parseable or not all mandatory arguments given
		return -1;
	}
	
	  // check for help request
	if (opts.getBoolVal("help")) {
		opts.coutUsage();
		return 0;
	}
	if (opts.getBoolVal("version")) {
		giveVersion();
		return 0;
	}

	  // get hydrophobicity scale table to use
	const double* scale = NULL;
	std::string scaleName = "none";
	switch (opts.getCharVal("scale")) {
		case 'K': scale = KyteDoolittle; scaleName="Kyte-Doolittle"; break;
		case 'G': scale = EngelmanGES; scaleName="GES/Engelman"; break;
		case 'E': scale = Eisenberg; scaleName="Eisenberg"; break;
		case 'H': scale = HoppWoodsNegated; scaleName="Hopp-Woods (negated)"; break;
		case 'C': scale = Cornette; scaleName="Cornette"; break;
		case 'J': scale = Janin; scaleName="Janin"; break;
		case 'R': scale = Rose; scaleName="Rose"; break;
		default:
			std::cerr <<"\n ERROR: given scale '" <<opts.getCharVal("scale")
						<<"' is not supported !\n\n";
			return -1;
	}
	
	double scaleAvg = 0.0;
	  // get mapping of amino acid to table position
	std::map< char, size_t > aa2pos;
	for (size_t i=0; i<aminoAcidNumber; i++) {
		aa2pos[aminoAcidCode[i]] = i;
		scaleAvg += scale[i];
	}
	  // get average scale
	scaleAvg /= (double)aminoAcidNumber;
	
		  // print used scale
	if (!opts.getBoolVal("s")) {
		std::cout <<"\n# used scale : " <<scaleName <<std::endl;
		if (opts.getBoolVal("v")) {
			std::cout <<"\n";
			double scaleAvgAbove = 0.0;
			size_t scaleAvgAboveCount = 0;
			double scaleAvgBelow = 0.0;
			for (size_t i=0; i<aminoAcidNumber; i++) {
				std::cout <<" " <<aminoAcidCode[i] <<" : " <<scale[i] <<"\n";
				if (scale[i] > scaleAvg) {
					scaleAvgAbove += scale[i];
					scaleAvgAboveCount++;
				} else {
					scaleAvgBelow += scale[i];
				}
			}
			scaleAvgAbove /= (double)scaleAvgAboveCount;
			scaleAvgBelow /= (double)(aminoAcidNumber - scaleAvgAboveCount);
			std::cout	<<"\n average of all : " <<scaleAvg 
						<<"\n average above  : " <<scaleAvgAbove 
						<<"\n average below  : " <<scaleAvgBelow <<"\n" <<std::endl;
		}
	}
	
	  // get amino acid sequence
	std::string aaSeq = opts.getStrVal("aa");
	if (aaSeq.size() == 0) {
		std::cerr <<"\n ERROR : no amino acid sequence given !\n\n";
		return -1;
	}
	  // check sequence alphabet
	for (size_t i=0; i<aaSeq.size(); i++) {
		if (aa2pos.find(aaSeq[i]) == aa2pos.end()) {
			std::cerr <<"\n ERROR : not supported character '" <<aaSeq[i] 
			        <<"' in amino acid sequence at pos. " 
					<<(i+1) <<" !\n\n";
			return -1;
		}
	}
	
	  // get and check the window size to use
	size_t winSize = 7;
	{
		int t= opts.getIntVal("size");
		if (t<1) {
			std::cerr <<"\n ERROR : the window size has to be at least 1 !\n\n";
			return -1;
		}
		if (t % 2 == 0) {
			std::cerr <<"\n ERROR : the window size has to be odd !\n\n";
			return -1;
		}
		if (t > (int)aaSeq.size()) {
			std::cerr <<"\n ERROR : the window size exeeds the sequence length !\n\n";
			return -1;
		}
		winSize = (size_t)t;
	}
	
	  // get minimal value to mark a window hydrophobic
	const double minHydro = opts.getDoubleVal("min");
	
	
	//////////////  CALCULATE AVERAGE VALUES  /////////////////////////////
	
	std::vector<double> hydro(aaSeq.size(),0.0);
	
	double curAvg = 0.0;
	size_t curPos = 0;
	size_t fillPos = 0;
	const size_t winSizeHalf = winSize / 2;
	 // generate values for the leading positions until a full window size is reached
	if (winSize > 1) {
		for (; curPos < winSizeHalf; curPos++) {
			curAvg += scale[aa2pos[aaSeq[curPos]]];
		}
		for (; curPos < winSize; curPos++ ) {
			curAvg += scale[aa2pos[aaSeq[curPos]]];
			hydro[ fillPos ] = ( curAvg / (double) (curPos+1) );
			fillPos++;
		}
	} else {
		curAvg += scale[aa2pos[aaSeq[curPos]]];
		curPos++;
		hydro[ fillPos ] = ( curAvg / (double) winSize );
		fillPos++;
	}
	 // handle full window sizes
	for (; curPos < aaSeq.size(); curPos++) {
		curAvg -= scale[aa2pos[aaSeq[curPos-winSize]]];
		curAvg += scale[aa2pos[aaSeq[curPos]]];
		hydro[ fillPos ] = ( curAvg / (double) winSize );
		fillPos++;
	}
	 // handle tailing decreasing window sizes
	if (winSize > 1) {
		for (;fillPos < aaSeq.size(); fillPos++) {
			curAvg -= scale[aa2pos[aaSeq[fillPos-winSizeHalf-1]]];
			hydro[ fillPos ] = ( curAvg / (double) (aaSeq.size()-fillPos+winSizeHalf) );
		}
	}
	
	if (!opts.getBoolVal("s")) {
		std::cout	<<"\n# averaged hydrophobicity values with window size "
					<<winSize <<" : \n" <<std::endl;
		std::cout <<hydro[0];
		for (size_t i=1; i<hydro.size(); i++) {
			std::cout <<" " <<hydro[i];
		}
		std::cout <<std::endl;
	}
	
	/////////////  DERIVE HP SEQUENCE FROM AVERAGE VALUES  /////////////////
	
	  // decide H or P based on the given minimal hydrophobicity value 
	std::string hpSeq = aaSeq;
	for (size_t i=0; i<aaSeq.size(); i++) {
		if (hydro[i] > minHydro) {
			hpSeq[i] = 'H';
		} else {
			hpSeq[i] = 'P';
		}
	}

	  // print HP sequence
	if (!opts.getBoolVal("s")) {
		std::cout	<<"\n# resulting HP sequence for min. hydrophobicity "
					<<minHydro <<" : \n\n";
	}
	std::cout	<<hpSeq
				<<std::endl;

	
	if (!opts.getBoolVal("s")) {
		std::cout	<<std::endl;
	}
	
	return 0;
}

									
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

void initAllowedArguments(	biu::OptionMap & allowedArgs,
							std::string &infoText ) 
{
	allowedArgs.push_back(biu::COption(	
							"aa", false, biu::COption::STRING, 
							"the amino acid sequence to get an HP sequence of"));
	allowedArgs.push_back(biu::COption(	
							"scale", true, biu::COption::CHAR, 
							"the hydrophobicity scale to use: (K)yte-Doolittle, (G)ES-Engelman, (E)isenberg, (H)opp-Woods, (C)ornette, (J)anin, or (R)ose"
							,"K"));
	allowedArgs.push_back(biu::COption(	
							"min", true, biu::COption::DOUBLE, 
							"minimal average value of a windows hydrophibicity to mark a position as hydrophobic"
							,"1.6"));
	allowedArgs.push_back(biu::COption(	
							"size", true, biu::COption::INT, 
							"the size of the window used to calculate the average for every position (has to be odd)"
							,"7"));
	allowedArgs.push_back(biu::COption(	
							"s", true, biu::COption::BOOL, 
							"silent output : only the resulting HP sequence"));
	allowedArgs.push_back(biu::COption(	
							"v", true, biu::COption::BOOL, 
							"verbose output"));
	allowedArgs.push_back(biu::COption(	
							"help", true, biu::COption::BOOL, 
							"program parameters and help"));
	allowedArgs.push_back(biu::COption(	
							"version", true, biu::COption::BOOL, 
							"version information of this program"));
	
	infoText = "\n"
			"Hydrophobicity scales available:\n"
			"\n"
			"Kyte-Doolittle scale: The Kyte-Doolittle scale is widely used for"
			" detecting hydrophobic regions in proteins. Regions with a"
			" positive value are hydrophobic. This scale can be used for"
			" identifying both surface-exposed regions as well as transmembrane"
			" regions, depending on the used window size. Short window sizes"
			" of 5-7 generally works well for predicting putative"
			" surface-exposed regions. Large window sizes of 19-21"
			" is well suited for finding transmembrane domains if the values"
			" calculated are above 1.6 [Kyte and Doolittle, 1982]."
			" These values should be used as a rule of thumb and deviations"
			" from the rule may occur.\n"
			"\n"
			"GES or Engelman scale: The Engelman hydrophobicity scale, also"
			" known as the GES-scale, is another scale which can be used for"
			" prediction of protein hydrophobicity [Engelman et al., 1986]."
			" As the Kyte-Doolittle scale, this scale is useful for"
			" prediction transmembrane regions in proteins.\n"
			"\n"
			"Eisenberg scale: The Eisenberg scale is a normalized consensus"
			" hydrophobicity scale, which share many features with the other"
			" hydrophobocity scales [Eisenberg et al., 1984].\n"
			"\n"
			"Cornette scale: Cornette et al. computed an optimal"
			" hydrophobicity scale based on 28 published scales"
			" [Cornette et al., 1987]. This optimized scale is also suitable"
			" for prediction of alpha-helices in proteins.\n"
			"\n"
			"Hopp-Woods scale: Hopp and Woods developed their hydrophobicity"
			" scale for identification of potential antigenic sites in"
			" proteins. This scale is basically a hydrophilic index where"
			" apolar residues have been assigned negative values. Antigenic"
			" sites are likely to be predicted when using a window size of"
			" 7 [Hopp and Woods, 1983].\n"
			" Note: Within this implementation all values are inverted to"
			" receive a hydrophobicity scale again!\n"
			"\n"
			"Rose scale: The hydrophobicity scale by Rose et al. is"
			" correlated to the average area of buried amino acids in"
			" globular proteins [Rose et al., 1985]. This results in a"
			" scale which is not showing the helices of a protein but rather"
			" the surface accessibility.\n"
			"\n"
			"Janin scale: This scale also tells about the accessible and"
			" buried amino acid residues of globular proteins [Janin, 1979].\n"
			"\n"
			"[source : http://www.clcbio.com]\n"
			"\n"
		;

}


////////////////////////////////////////////////////////////////////////////

//	    (1)   (2)   (3)   (4)   (5)   (6)   (7)
//	A   1.80 -0.50  0.20  0.62  0.74  0.30  1.60
//	C   2.50 -1.00  4.10  0.29  0.91  0.90  2.00
//	D  -3.50  3.00 -3.10 -0.90  0.62 -0.60 -9.20
//	E  -3.50  3.00 -1.80 -0.74  0.62 -0.70 -8.20
//	F   2.80 -2.50  4.40  1.19  0.88  0.50  3.70
//	G  -0.40  0.00  0.00  0.48  0.72  0.30  1.00
//	H  -3.20 -0.50  0.50 -0.40  0.78 -0.10 -3.00
//	I   4.50 -1.80  4.80  1.38  0.88  0.70  3.10
//	K  -3.90  3.00 -3.10 -1.50  0.52 -1.80 -8.80
//	L   3.80 -1.80  5.70  1.06  0.85  0.50  2.80
//	M   1.90 -1.30  4.20  0.64  0.85  0.40  3.40
//	N  -3.50  0.20 -0.50 -0.78  0.63 -0.50 -4.80
//	P  -1.60  0.00 -2.20  0.12  0.64 -0.30 -0.20
//	Q  -3.50  0.20 -2.80 -0.85  0.62 -0.70 -4.10
//	R  -4.50  3.00  1.40 -2.53  0.64 -1.40 -12.3
//	S  -0.80  0.30 -0.50 -0.18  0.66 -0.10  0.60
//	T  -0.70 -0.40 -1.90 -0.05  0.70 -0.20  1.20
//	V   4.20 -1.50  4.70  1.08  0.86  0.60  2.60
//	W  -0.90 -3.40  1.00  0.81  0.85  0.30  1.90
//	Y  -1.30 -2.30  3.20  0.26  0.76 -0.40 -0.70
//	
//	(1) Kyte-Doolittle scale: The Kyte-Doolittle scale is widely used for detecting hydrophobic regions in proteins. Regions with a positive value are hydrophobic. This scale can be used for
//	identifying both surface-exposed regions as well as transmembrane regions, depending on the used window size. Short window sizes of 5-7 generally works well for predicting putative
//	surface-exposed regions. Large window sizes of 19-21 is well suited for finding transmembrane domains if the values calculated are above 1.6 [Kyte and Doolittle, 1982]. These values
//	should be used as a rule of thumb and deviations from the rule may occur.
//	
//	(7) Engelman scale: The Engelman hydrophobicity scale, also known as the GES-scale, is another scale which can be used for prediction of protein hydrophobicity [Engelman et al., 1986].
//	As the Kyte-Doolittle scale, this scale is useful for prediction transmembrane regions in proteins.
//	
//	(4) Eisenberg scale: The Eisenberg scale is a normalized consensus hydrophobicity scale, which share many features with the other hydrophobocity scales [Eisenberg et al., 1984].
//	
//	(2) Hopp-Woods scale: Hopp and Woods developed their hydrophobicity scale for identification of potential antigenic sites in proteins. This scale is basically a hydrophilic index where apolar
//	residues have been assigned negative values. Antigenic sites are likely to be predicted when using a window size of 7 [Hopp and Woods, 1983].
//	
//	(3) Cornette scale: Cornette et al. computed an optimal hydrophobicity scale based on 28 published scales [Cornette et al., 1987]. This optimized scale is also suitable for prediction of
//	alpha-helices in proteins.
//	
//	(5) Rose scale: The hydrophobicity scale by Rose et al. is correlated to the average area of buried amino acids in globular proteins [Rose et al., 1985]. This results in a scale which is not
//	showing the helices of a protein but rather the surface accessibility.
//	
//	(6) Janin scale: This scale also tells about the accessible and buried amino acid residues of globular proteins [Janin, 1979].

//	source : http://www.clcbio.com
