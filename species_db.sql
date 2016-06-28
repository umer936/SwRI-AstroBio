-- Species,+/-,H,He,C,N,O,F,P,S,Cl,Na,Mg,Al,Si,K,Ca,Ti,Fe,MolMass,n,f,T3,[K],P3,[kPa]

-- CREATE DATABASE astrobio;

DROP TABLE species_all;
DROP TABLE species;

CREATE TABLE species_all (id MEDIUMINT UNSIGNED NOT NULL AUTO_INCREMENT,
                      name VARCHAR (20),
                      charge TINYINT default 0, 
                      hydrogen TINYINT UNSIGNED default 0,
                      helium TINYINT UNSIGNED default 0,
                      carbon TINYINT UNSIGNED default 0,
                      nitrogen TINYINT UNSIGNED default 0,
                      oxygen TINYINT UNSIGNED default 0,
                      fluorine TINYINT UNSIGNED default 0,
                      phosphorus TINYINT UNSIGNED default 0,
                      sulfur TINYINT UNSIGNED default 0,
                      chlorine TINYINT UNSIGNED default 0,
                      sodium TINYINT UNSIGNED default 0,
                      magnesium TINYINT UNSIGNED default 0,
                      aluminum TINYINT UNSIGNED default 0,
                      silicon TINYINT UNSIGNED default 0,
                      potassium TINYINT UNSIGNED default 0,
                      calcium TINYINT UNSIGNED default 0,
                      titanium TINYINT UNSIGNED default 0,
                      iron TINYINT UNSIGNED default 0,
                      molmass float,
                      n TINYINT,
                      f TINYINT,
                      temperature float,
                      pressure float,
                      PRIMARY KEY (id)
                     );
CREATE TABLE species (id MEDIUMINT UNSIGNED NOT NULL AUTO_INCREMENT,
                      name VARCHAR (20),
                      charge TINYINT default 0, 
                      hydrogen TINYINT UNSIGNED default 0,
                      helium TINYINT UNSIGNED default 0,
                      carbon TINYINT UNSIGNED default 0,
                      nitrogen TINYINT UNSIGNED default 0,
                      oxygen TINYINT UNSIGNED default 0,
                      fluorine TINYINT UNSIGNED default 0,
                      phosphorus TINYINT UNSIGNED default 0,
                      sulfur TINYINT UNSIGNED default 0,
                      chlorine TINYINT UNSIGNED default 0,
                      sodium TINYINT UNSIGNED default 0,
                      magnesium TINYINT UNSIGNED default 0,
                      aluminum TINYINT UNSIGNED default 0,
                      silicon TINYINT UNSIGNED default 0,
                      potassium TINYINT UNSIGNED default 0,
                      calcium TINYINT UNSIGNED default 0,
                      titanium TINYINT UNSIGNED default 0,
                      iron TINYINT UNSIGNED default 0,
                      molmass float,
                      n TINYINT,
                      f TINYINT,
                      temperature float,
                      pressure float,
                      UNIQUE (charge, hydrogen, helium, carbon, nitrogen, fluorine, phosphorus, sulfur, chlorine, sodium, magnesium, aluminum, silicon, potassium, calcium, titanium),
                      PRIMARY KEY (id)
                     );

LOAD DATA INFILE '/web/phidrates/AstroBio/SpeciesMNew.txt' IGNORE INTO TABLE species_all FIELDS TERMINATED BY ',' IGNORE 1 LINES 
    (name, charge, hydrogen, helium, carbon, nitrogen, oxygen, fluorine, phosphorus, sulfur, chlorine, sodium, magnesium, aluminum, silicon, potassium, calcium, titanium, iron, molmass, n, f, temperature, pressure)
SET id=NULL;
LOAD DATA INFILE '/web/phidrates/AstroBio/SpeciesMNew.txt' IGNORE INTO TABLE species FIELDS TERMINATED BY ',' IGNORE 1 LINES 
    (name, charge, hydrogen, helium, carbon, nitrogen, oxygen, fluorine, phosphorus, sulfur, chlorine, sodium, magnesium, aluminum, silicon, potassium, calcium, titanium, iron, molmass, n, f, temperature, pressure)
SET id=NULL;

SELECT name, charge, hydrogen, helium, carbon, nitrogen, fluorine, phosphorus, sulfur, chlorine, sodium, magnesium, aluminum, silicon, potassium, calcium, titanium, iron, molmass, n, f, temperature, pressure FROM species INTO OUTFILE '/tmp/SpeciesNoDupes.txt';
SELECT name, charge, hydrogen, helium, carbon, nitrogen, fluorine, phosphorus, sulfur, chlorine, sodium, magnesium, aluminum, silicon, potassium, calcium, titanium, iron, molmass, n, f, temperature, pressure FROM species_all INTO OUTFILE '/tmp/SpeciesAll.txt';
