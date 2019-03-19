tsir_track_list <- list(
	finkenstadt2000time=data.frame(
		author=c("Finkenstadt", "Grenfell"),
		year=2000,
		method="Regression",
		disease="Measles"
	),
	bjornstad2002dynamics=data.frame(
		author=c("Bjornstad", "Finkenstadt", "Grenfell"),
		year=2002,
		method="Regression",
		disease="Measles"
	),
	grenfell2002dynamics=data.frame(
		author=c("Grenfell", "Bjornstad", "Finkenstadt"),
		year=2002,
		method="Regression",
		disease="Measles"
	),
	finkenstadt2002stochastic=data.frame(
		author=c("Finkenstadt", "Bjornstad", "Grenfell"),
		year=2002,
		method="Regression",
		disease="Measles"
	),
	clark2004population=data.frame(
		author=c("Clark", "Bjornstad"),
		year=2004,
		method="MCMC",
		disease="Measles"
	),
	morton2005discrete=data.frame(
		author=c("Morton", "Finkenstadt"),
		year=2005,
		method="MCMC",
		disease="Measles"
	),
	ferrari2008dynamics=data.frame(
		author=c("Ferrari", "Grais", "Bharti", "Conlan", "Bjornstad", "Wolfson", "Guerin", "Djibo", "Grenfell"),
		year=2008,
		method="MCMC",
		disease="Measles"
	),
	koelle2004disentangling=data.frame(
		author=c("Koelle", "Pascual"),
		year=2004,
		method="Regression",
		disease="Cholera"
	),
	pascual2007shifting=data.frame(
		author=c("Pascual", "Cazelles", "Bouma", "Chaves", "Koelle"),
		year=2007,
		method="Regression",
		disease="Malaria"
	),
	metcalf2009seasonality=data.frame(
		author=c("Metcalf", "Grenfell"),
		year=2009,
		method="Regression",
		disease=c("Measles", "Pertussis", "mumps", "diphtheria", "Varicella", "Scarlet fever")
	),
	pascual2008predicting=data.frame(
		author=c("Pascual", "Chaves", "Cash", "Rodo", "Yunus"),
		year=2008,
		method="Regression",
		disease="Cholera"
	),
	takahash2016hand=data.frame(
		author=c("Takahashi", "Liao", "Van Boeckel", "Xing", "Sun", "Hsiao", "Metcalf", "Chang", 
				 "Liu", "Zhang", "Wu", "Cowling", "Leung", "Farrar", "Rogier van Doorn", "Grenfell", "Yu"),
		year=2016,
		method="Regression",
		disease="HFMD"
	),
	metcalf2010rubella=data.frame(
		author=c("Metcalf", "Munayco", "Chowell", "Grenfell", "Bjornstad"),
		year=2010,
		method="Regression",
		disease="Rubella"
	),
	ferrari2010rural=data.frame(
		author=c("Ferrari", "Djibo", "Grais", "Bharti", "Grenfell", "Bjornstad"),
		year=2010,
		method="Regression",
		disease="Measles"
	),
	perkins2015estimating=data.frame(
		author=c("Perkins", "Metcalf", "Grenfell", "Tatem"),
		year=2015,
		method="Regression",
		disease="Chikungunya virus"
	),
	metcalf2011epidemiology=data.frame(
		author=c("Metcalf", "Bjornstad", "Ferrari", "Klepac", "Bharti", "Lopez-Gatell", "Grenfell"),
		year=2011,
		method="MCMC",
		disease="Rubella"
	),
	dalziel2016persistent=data.frame(
		author=c("Dalziel", "Bjornstad", "van Panhuis", "Burke", "Metcalf", "Grenfell"),
		year=2016,
		method="Regression",
		disease="Measles"
	),
	metcalf2013implications=data.frame(
		author=c("Metcalf", "Cohen", "Lessler", "McAnerney", "Ntshoe", "Puren", "Klepac", "Tatem",
				 "Grenfell", "Bjorstad"),
		year=2013,
		method="Regression",
		disease="Rubella"
	),
	jackson2014effects=data.frame(
		author=c("Jackson", "Mangtani", "Fine", "Vynnycky"),
		year=2014,
		method="Regression",
		disease="Varicella"
	),
	tian2017anthropogenically=data.frame(
		author=c("Tian", "Yu", "Bjornstad", "Cazelles", "Yang", "Tan", "Huang", "Cui", "Dong", "Ma",
				 "Ma", "Zhou", "Laine", "Wu", "Zhang", "Wang", "Yang", "Stenseth", "Xu"),
		year=2017,
		method="MCMC",
		disease="Hemorrhagic fever"
	),
	joh2013dynamics=data.frame(
		author=c("Joh", "Hoekstra", "Barzilay", "Bowen", "Mintz", "Weiss", "Weitz"),
		year=2013,
		method="MCMC",
		disease="Shigellosis"
	),
	mantilla2009decreasing=data.frame(
		author=c("Mantilla-Beniers", "Bjornstad", "Grenfell", "Rohani"),
		year=2009,
		method="Regression",
		disease="Measles"
	),
	caudron2015predictability=data.frame(
		author=c("Caudron", "Mahmud", "Metcalf", "Gottfreosson", "Viboud", "Cliff", "Grenfell"),
		year=2015,
		method="Regression",
		disease="Measles"
	),
	van2016hand=data.frame(
		author=c("Van Boeckel", "Takahashi", "Liao", "Xing", "Lai", "Hsiao", "Liu", "Zheng", "Chang",
				 "Yuan", "Metcalf", "Yu", "Grenfell"),
		year=2016,
		method="Regression",
		disease="HFMD"
	),
	mahmud2017comparative=data.frame(
		author=c("Mahmud", "Metcalf", "Grenfell"),
		year=2017,
		method="Regression",
		disease=c("Measles", "Pertussis", "mumps", "diphtheria", "Varicella", "Scarlet fever")
	),
	kraemer2018inferences=data.frame(
		author=c("Kraemer", "Bisanzio", "Reiner", "Zakar", "Hawkins", "Freifeld", "Smith", "Hay",
				 "Brownstein", "Perkins"),
		year=2018,
		method="Regression",
		disease="Dengue"
	),
	takahashi2018epidemic=data.frame(
		author=c("Takahashi", "Metcalf", "Arima", "Fujimoto", "Shimizu", "Rogier van Doorn", "Le Van",
				 "Chan", "Farrar", "Oishi", "Grenfell"),
		year=2018,
		method="Regression",
		disease="HFMD"
	),
	du2017estimating=data.frame(
		author=c("Du", "Zhang", "Zhang", "Yu", "Hao"),
		year=2017,
		method="Regression",
		disease="HFMD"
	),
	harris2018climate=data.frame(
		author=c("Harris", "Caldwell", "Mordecai"),
		year=2018,
		method="Regression",
		disease="Zika"
	),
	baker2018dynamic=data.frame(
		author=c("Baker", "Mahmud", "Metcalf"),
		year=2018,
		method="Regression",
		disease="Varicella"
	),
	mahmud2017drivers=data.frame(
		author=c("Mahmud", "Alam", "Metcalf"),
		year=2017,
		method="Regression",
		disease="Measles"
	),
	reich2013interactions=data.frame(
		author=c("Reich", "Shrestha", "King", "Rohani", "Lessler", "Kalayanarooj", "Yoon", "Gibbons", 
				 "Burke", "Cummings"),
		year=2013,
		method="Regression",
		disease="Dengue"
	),
	becker2017tsir=data.frame(
		author=c("Becker", "Grenfell"),
		year=2017,
		method="Regression",
		disease="Measles"
	)
)
