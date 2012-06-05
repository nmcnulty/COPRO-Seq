package genomecodes;

# Eventually need to automate updates to this module such that a 3-column table can be used as 
# input, where column 1 is accession #, 2 is abbreviation, 3 is species name

sub declare_inhouseacc_for_genomeabbrev {
	my %inhouseacc_for_genomeabbrev = (
		"BACCAC"	=> "NZ_AAVM00000000",
		"BACOVA" 	=> "NZ_AAXF00000000",
		"BACTHE7330"	=> "NC_Bthetaiotaomicron7330",
		"BACTHE" 	=> "NC_004663",
		"BACUNI" 	=> "NZ_AAYH00000000",
		"BACVUL" 	=> "NC_009614",
		"BACWH2" 	=> "NC_BWH2",
		"BACWH2-CLOSED"	=> "NC_BacWH2_v1",
		"BIFANI" 	=> "NZ_BanimalisDN1730010",
		"BLAHYD"	=> "NZ_ACBZ00000000",
		"BRYFOR"	=> "NZ_ACCL00000000",
		"CLOSCI" 	=> "NZ_ABFY00000000",
		"CLOSPI" 	=> "NZ_ABIK00000000",
		"CLOSYM"	=> "NC_Csymbiosum",
		"COLAER" 	=> "NZ_AAVN00000000",
		"DESPIG"	=> "NZ_ABXU00000000",
		"DESPIGGOR1"	=> "NC_DPigerGOR1",
		"DORLON" 	=> "NZ_AAXB00000000",
		"DORFOR"	=> "NZ_AAXA00000000",
		"ESCCOL"	=> "NC_000913",
		"EUBHAL"	=> "NZ_ACEP00000000",
		"EUBREC" 	=> "NC_012781",
		"FAEPRA" 	=> "NZ_ABED00000000",
		"LACDEL130"	=> "NC_LB130",
		"LACDEL182" 	=> "NC_LB182",
		"LACLAC066"	=> "NC_LACLAC",
		"METSMI"	=> "NC_009515",
		"PARDIS" 	=> "NC_009615",
		"ROSINT"	=> "NZ_ABYJ00000000",
		"RUMOBE" 	=> "NZ_AAVO00000000",
		"RUMTOR" 	=> "NZ_AAVP00000000",
		"STRTHE147"	=> "NC_SthermophilusDN001147",			# NOT the Activia strain
		"STRTHE171" => "NZ_SthermophilusDN001171"			# The Activia strain
	);
	return %inhouseacc_for_genomeabbrev;
}

sub declare_speciesname_for_inhouseacc {
	my %speciesname_for_inhouseacc = (
		"NZ_AAVM00000000"			=> "B_caccae_ATCC_43185",
		"NZ_ACCL00000000"			=> "B_formatexigens_DSM_14469",
		"NZ_ACBZ00000000"			=> "B_hydrogenotrophica_DSM_10507",
		"NZ_AAXF00000000"			=> "B_ovatus_ATCC_8483",
		"NC_004663"					=> "B_thetaiotaomicron_VPI-5482",
		"NC_Bthetaiotaomicron7330"		=> "B_thetaiotaomicron_7330",
		"NZ_AAYH00000000"			=> "B_uniformis_ATCC_8492",
		"NC_009614"					=> "B_vulgatus_ATCC_8482",
		"NC_BWH2"					=> "B_WH2",
		"NC_BacWH2_v1"			=> "B_WH2_closed",
		"NZ_BanimalisDN1730010"		=> "B_animalis_subsp_lactis_CNCM_I-2494",
		"NZ_ABFY00000000"			=> "C_scindens_ATCC_35704",
		"NZ_ABIK00000000"			=> "C_spiroforme_DSM_1552",
		"NC_Csymbiosum"				=> "C_symbiosum",
		"NZ_AAVN00000000"			=> "C_aerofaciens_ATCC_25986",
		"NZ_AAXA00000000"			=> "D_formicigenerans",
		"NZ_AAXB00000000"			=> "D_longicatena_DSM_13814",
		"NZ_ABXU00000000"			=> "D_piger_ATCC_29098",
		"NC_DpigerGOR1"			=> "D_piger_GOR1",
		"NC_000913"			=> "E_coli_K-12_substr_MG1655",
		"NZ_ACEP00000000"			=> "E_hallii",
		"NC_012781"					=> "E_rectale_ATCC_33656",
		"NZ_ABED00000000"			=> "F_prausnitzii_M21-2",
		"NC_LB130"	=> "L_delbrueckii_subsp_bulgaricus_CNCM_I-1632",
		"NC_LB182"	=> "L_delbrueckii_subsp_bulgaricus_CNCM_I-1519",
		"NC_LACLAC"		=> "L_lactis_subsp_cremoris_CNCM_I-1631",
		"NC_009515"					=> "M_smithii_ATCC_35061",
		"NC_009615"					=> "P_distasonis_ATCC_8503",
		"NZ_ABYJ00000000"			=> "R_intestinalis_L1",
		"NZ_AAVO00000000"			=> "R_obeum_ATCC_29174",
		"NZ_AAVP00000000"			=> "R_torques_ATCC_27756",
		"NC_SthermophilusDN001147"	=> "S_thermophilus_DN-001147",
		"NZ_SthermophilusDN001171"	=> "S_thermophilus_CNCM_I-1630"
	);
	
	return %speciesname_for_inhouseacc;
}

# NCBI codes below are not being maintained until I can figure out how to get complete genomes pulled down from GenBank

sub declare_genbankacc_for_genomeabbrev {
	my %genbankacc_for_genomeabbrev = (
		"BACCAC"		=> "NZ_AAVM00000000",
		"BACOVA"		=> "NZ_AAXF00000000",
		"BACTHE"		=> "NC_004663",
		"BACUNI"		=> "NZ_AAYH00000000",
		"BACVUL"		=> "NC_009614",
		# No GenBank entry for BACWH2 yet...
		# No GenBank entry for B. animalis DN 173-010 yet...
		"CLOSCI"		=> "NZ_ABFY00000000",
		"CLOSPI"		=> "NZ_ABIK00000000",
		"COLAER"		=> "NZ_AAVN00000000",
		"DORLON"		=> "NZ_AAXB00000000",
		"EUBREC"		=> "NC_012781",
		"FAEPRA"		=> "NZ_ABED00000000",
		# No GenBank entry for L. delbrueckii bulgaricus DN 100-130 yet...
		# No GenBank entry for L. delbrueckii bulgaricus DN 100-182 yet...
		# No GenBank entry for L. lactis DN 030-066 yet...
		"PARDIS"		=> "NC_009615",
		"RUMOBE"		=> "NZ_AAVO00000000",
		"RUMTOR"		=> "NZ_AAVP00000000"
		# No GenBank entry for S. thermophilus DN 001-147 or 001-171 yet...
	);
	return %genbankacc_for_genomeabbrev;
}

sub declare_speciesname_for_genbankacc {
	my %speciesname_for_genbankacc = (
		"NZ_AAVM00000000"		=> "B_caccae",
		"NZ_AAXF00000000"		=> "B_ovatus",
		"NC_004663"				=> "B_thetaiotaomicron",
		"NZ_AAYH00000000"		=> "B_uniformis",
		"NC_009614"				=> "B_vulgatus",
		# No GenBank entry for BACWH2 yet...
		# No GenBank entry for B. animalis DN 173-010 yet...
		"NZ_ABFY00000000"		=> "C_scindens",
		"NZ_ABIK00000000"		=> "C_spiroforme",
		"NZ_AAVN00000000"		=> "C_aerofaciens",
		"NZ_AAXB00000000"		=> "D_longicatena",
		"NC_012781"				=> "E_rectale",
		"NZ_ABED00000000"		=> "F_prausnitzii M21-2",
		# No GenBank entry for L. delbrueckii bulgaricus DN 100-130 yet...
		# No GenBank entry for L. delbrueckii bulgaricus DN 100-182 yet...
		# No GenBank entry for L. lactis DN 030-066 yet...
		"NC_009615"				=> "P_distasonis",
		"NZ_AAVO00000000"		=> "R_obeum",
		"NZ_AAVP00000000"		=> "R_torques"
		# No GenBank entry for S. thermophilus DN 001-147 or 001-171 yet...
	);
	return %speciesname_for_genbankacc;
}

1;
