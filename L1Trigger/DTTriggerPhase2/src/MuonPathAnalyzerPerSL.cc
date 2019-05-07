#include "L1Trigger/DTTriggerPhase2/interface/MuonPathAnalyzerPerSL.h"

using namespace edm;
using namespace std;



// ============================================================================
// Constructors and destructor
// ============================================================================
MuonPathAnalyzerPerSL::MuonPathAnalyzerPerSL(const ParameterSet& pset) :
  MuonPathAnalyzer(pset),
  bxTolerance(30),
  minQuality(LOWQGHOST),
  chiSquareThreshold(50)
{
  // Obtention of parameters
  debug         = pset.getUntrackedParameter<Bool_t>("debug");
  if (debug) cout <<"MuonPathAnalyzerPerSL: constructor" << endl;

  //z
  int rawId;
  z_filename = pset.getUntrackedParameter<std::string>("z_filename");
  std::ifstream ifin2(z_filename.c_str());
  double z;
  while (ifin2.good()){
    ifin2 >> rawId >> z;
    zinfo[rawId]=z;
  }
  
  //shift
  shift_filename = pset.getUntrackedParameter<std::string>("shift_filename");
  std::ifstream ifin3(shift_filename.c_str());
  double shift;
  while (ifin3.good()){
    ifin3 >> rawId >> shift;
    shiftinfo[rawId]=shift;
  }
  
  chosen_sl = pset.getUntrackedParameter<int>("trigger_with_sl");
    
  if(chosen_sl!=1 && chosen_sl!=3 && chosen_sl!=4){
    std::cout<<"chosen sl must be 1,3 or 4(both superlayers)"<<std::endl;
    assert(chosen_sl!=1 && chosen_sl!=3 && chosen_sl!=4); //4 means run using the two superlayers
  }
  
}


MuonPathAnalyzerPerSL::~MuonPathAnalyzerPerSL() {
  if (debug) cout <<"MuonPathAnalyzerPerSL: destructor" << endl;
}


// ============================================================================
// Main methods (initialise, run, finish)
// ============================================================================
void MuonPathAnalyzerPerSL::initialise(const edm::EventSetup& iEventSetup) {
  if(debug) cout << "MuonPathAnalyzerPerSL::initialiase" << endl;

  if(debug) std::cout<<"getting DT geometry"<<std::endl;
  iEventSetup.get<MuonGeometryRecord>().get(dtGeo);//1103

  
}


void MuonPathAnalyzerPerSL::run(edm::Event& iEvent, const edm::EventSetup& iEventSetup, std::vector<MuonPath*> &muonpaths, std::vector<metaPrimitive> &metaPrimitives) {
  if (debug) cout <<"MuonPathAnalyzerPerSL: run" << endl;
  
  // fit per SL (need to allow for multiple outputs for a single mpath) 
  for(auto muonpath = muonpaths.begin();muonpath!=muonpaths.end();++muonpath) {
    
    analyze(*muonpath, metaPrimitives);   
    
  }
  
  // Once we have the fit per SL, try to merge 
  //  associate(mpaths);  

}

void MuonPathAnalyzerPerSL::run(edm::Event& iEvent, const edm::EventSetup& iEventSetup, std::vector<MuonPath*> &muonpaths, std::vector<MuonPath*> &outMPath) {
  if (debug) cout <<"MuonPathAnalyzerPerSL: run" << endl;
  
  // fit per SL (need to allow for multiple outputs for a single mpath) 
  for(auto muonpath = muonpaths.begin();muonpath!=muonpaths.end();++muonpath) {
    
    analyze(*muonpath, outMPath);   
    
  }
  
  // Once we have the fit per SL, try to merge 
  //  associate(mpaths);  

}

void MuonPathAnalyzerPerSL::finish() {
  if (debug) cout <<"MuonPathAnalyzerPerSL: finish" << endl;
};


const int MuonPathAnalyzerPerSL::LAYER_ARRANGEMENTS[MAX_VERT_ARRANG][3] = {
    {0, 1, 2}, {1, 2, 3},                       // Grupos consecutivos
    {0, 1, 3}, {0, 2, 3}                        // Grupos salteados
};

// ============================================================================
// Private (dedicated) methods 
// ============================================================================

void MuonPathAnalyzerPerSL::setCellLayout(const int layout[4]) {
    memcpy(cellLayout, layout, 4 * sizeof(int));
    //celllayout[0]=layout[0];
    //celllayout[1]=layout[1];
    //celllayout[2]=layout[2];
    //celllayout[3]=layout[3];
    
    buildLateralities();
}

void MuonPathAnalyzerPerSL::buildLateralities(void) {
  
  LATERAL_CASES (*validCase)[4], sideComb[4];
  
  totalNumValLateralities = 0;
  // We generate all the possible laterality combinations for the group
  // of cells that belong to a given mpath.
  /* Generamos todas las posibles combinaciones de lateralidad para el grupo
     de celdas que forman parte del analizador */
  for(int lowLay = LEFT; lowLay <= RIGHT; lowLay++) {
    for(int midLowLay = LEFT; midLowLay <= RIGHT; midLowLay++) {
      for(int midHigLay = LEFT; midHigLay <= RIGHT; midHigLay++) {
	for(int higLay = LEFT; higLay <= RIGHT; higLay++) {
	  
	  sideComb[0] = static_cast<LATERAL_CASES>(lowLay);
	  sideComb[1] = static_cast<LATERAL_CASES>(midLowLay);
	  sideComb[2] = static_cast<LATERAL_CASES>(midHigLay);
	  sideComb[3] = static_cast<LATERAL_CASES>(higLay);
	  
	  /* Si una combinación de lateralidades es válida, la almacenamos */
	  if (isStraightPath(sideComb)) {
	    validCase = lateralities + totalNumValLateralities;
	    memcpy(validCase, sideComb, 4 * sizeof(LATERAL_CASES));
	    
	    latQuality[totalNumValLateralities].valid            = false;
	    latQuality[totalNumValLateralities].bxValue          = 0;
	    latQuality[totalNumValLateralities].quality          = NOPATH;
	    latQuality[totalNumValLateralities].invalidateHitIdx = -1;
	    
	    totalNumValLateralities++;
	  }
	}
      }
    }
  }
  
}
bool MuonPathAnalyzerPerSL::isStraightPath(LATERAL_CASES sideComb[4]) {

    return true; //trying with all lateralities to be confirmed

    int i, ajustedLayout[4], pairDiff[3], desfase[3];

    /* Sumamos el valor de lateralidad (LEFT = 0, RIGHT = 1) al desfase
       horizontal (respecto de la celda base) para cada celda en cuestion */
    for(i = 0; i <= 3; i++) ajustedLayout[i] = cellLayout[i] + sideComb[i];
    /* Variación del desfase por pares de celdas consecutivas */
    for(i = 0; i <= 2; i++) pairDiff[i] = ajustedLayout[i+1] - ajustedLayout[i];
    /* Variación de los desfases entre todas las combinaciones de pares */
    for(i = 0; i <= 1; i++) desfase[i] = abs(pairDiff[i+1] - pairDiff[i]);
    desfase[2] = abs(pairDiff[2] - pairDiff[0]);
    /* Si algún desfase es mayor de 2 entonces la trayectoria no es recta */
    bool resultado = (desfase[0] > 1 or desfase[1] > 1 or desfase[2] > 1);

    return ( !resultado );
}

void MuonPathAnalyzerPerSL::evaluatePathQuality(MuonPath *mPath) {

    int totalHighQ = 0, totalLowQ = 0;
  
    if(debug) std::cout<<"DTp2:evaluatePathQuality \t\t\t\t\t En evaluatePathQuality Evaluando PathQ. Celda base: "<<mPath->getBaseChannelId()<<std::endl;
    if(debug) std::cout<<"DTp2:evaluatePathQuality \t\t\t\t\t Total lateralidades: "<<totalNumValLateralities<<std::endl;

    // Por defecto.
    mPath->setQuality(NOPATH);

    /* Ensayamos los diferentes grupos de lateralidad válidos que constituyen
       las posibles trayectorias del muón por el grupo de 4 celdas.
       Posiblemente esto se tenga que optimizar de manera que, si en cuanto se
       encuentre una traza 'HIGHQ' ya no se continue evaluando mas combinaciones
       de lateralidad, pero hay que tener en cuenta los fantasmas (rectas
       paralelas) en de alta calidad que se pueden dar en los extremos del BTI.
       Posiblemente en la FPGA, si esto se paraleliza, no sea necesaria tal
       optimización */
    for (int latIdx = 0; latIdx < totalNumValLateralities; latIdx++) {
	if(debug) std::cout<<"DTp2:evaluatePathQuality \t\t\t\t\t Analizando combinacion de lateralidad: "
			   <<lateralities[latIdx][0]<<" "
			   <<lateralities[latIdx][1]<<" "	 
			   <<lateralities[latIdx][2]<<" "
			   <<lateralities[latIdx][3]<<std::endl;
	  
	evaluateLateralQuality(latIdx, mPath, &(latQuality[latIdx]));

	if (latQuality[latIdx].quality == HIGHQ) {
	    totalHighQ++;
	    if(debug) std::cout<<"DTp2:evaluatePathQuality \t\t\t\t\t\t Lateralidad HIGHQ"<<std::endl;
	}
	if (latQuality[latIdx].quality == LOWQ) {
	    totalLowQ++;
	    if(debug) std::cout<<"DTp2:evaluatePathQuality \t\t\t\t\t\t Lateralidad LOWQ"<<std::endl;
	}
    }
    /*
     * Establecimiento de la calidad.
     */
    if (totalHighQ == 1) {
	mPath->setQuality(HIGHQ);
    }
    else if (totalHighQ > 1) {
	mPath->setQuality(HIGHQGHOST);
    }
    else if (totalLowQ == 1) {
	mPath->setQuality(LOWQ);
    }
    else if (totalLowQ > 1) {
	mPath->setQuality(LOWQGHOST);
    }
}

void MuonPathAnalyzerPerSL::evaluateLateralQuality(int latIdx, MuonPath *mPath, LATQ_TYPE *latQuality){

    int layerGroup[3];
    LATERAL_CASES sideComb[3];
    PARTIAL_LATQ_TYPE latQResult[4] = {
	{false, 0}, {false, 0}, {false, 0}, {false, 0}
    };

    // Default values.
    latQuality->valid            = false;
    latQuality->bxValue          = 0;
    latQuality->quality          = NOPATH;
    latQuality->invalidateHitIdx = -1;

    /* En el caso que, para una combinación de lateralidad dada, las 2
       combinaciones consecutivas de 3 capas ({0, 1, 2}, {1, 2, 3}) fueran
       traza válida, habríamos encontrado una traza correcta de alta calidad,
       por lo que sería innecesario comprobar las otras 2 combinaciones
       restantes.
       Ahora bien, para reproducir el comportamiento paralelo de la FPGA en el
       que el análisis se va a evaluar simultáneamente en todas ellas,
       construimos un código que analiza las 4 combinaciones, junto con una
       lógica adicional para discriminar la calidad final de la traza */
    for (int i = 0; i <= 3 ; i++) {
	memcpy(layerGroup, LAYER_ARRANGEMENTS[i], 3 * sizeof(int));
	
	// Seleccionamos la combinación de lateralidad para cada celda.
	for (int j = 0; j < 3; j++)
	  sideComb[j] = lateralities[latIdx][ layerGroup[j] ];
	
	validate(sideComb, layerGroup, mPath, &(latQResult[i]));
    }
    /*
      Imponemos la condición, para una combinación de lateralidad completa, que
      todas las lateralidades parciales válidas arrojen el mismo valor de BX
      (dentro de un margen) para así dar una traza consistente.
      En caso contrario esa combinación se descarta.
    */
    if ( !sameBXValue(latQResult) ) {
	// Se guardan en los default values inciales.
	if(debug) std::cout<<"DTp2:evaluateLateralQuality \t\t\t\t\t Lateralidad DESCARTADA. Tolerancia de BX excedida"<<std::endl;
	return;
    }

    // Dos trazas complementarias válidas => Traza de muón completa.
    if ((latQResult[0].latQValid && latQResult[1].latQValid) or
	(latQResult[0].latQValid && latQResult[2].latQValid) or
	(latQResult[0].latQValid && latQResult[3].latQValid) or
	(latQResult[1].latQValid && latQResult[2].latQValid) or
	(latQResult[1].latQValid && latQResult[3].latQValid) or
	(latQResult[2].latQValid && latQResult[3].latQValid))
	{
	    latQuality->valid   = true;
	    //     latQuality->bxValue = latQResult[0].bxValue;
	    /*
	     * Se hace necesario el contador de casos "numValid", en vez de promediar
	     * los 4 valores dividiendo entre 4, puesto que los casos de combinaciones
	     * de 4 hits buenos que se ajusten a una combinación como por ejemplo:
	     * L/R/L/L, dan lugar a que en los subsegmentos 0, y 1 (consecutivos) se
	     * pueda aplicar mean-timer, mientras que en el segmento 3 (en el ejemplo
	     * capas: 0,2,3, y combinación L/L/L) no se podría aplicar, dando un
	     * valor parcial de BX = 0.
	     */
	    int sumBX = 0, numValid = 0;
	    for (int i = 0; i <= 3; i++) {
		if (latQResult[i].latQValid) {
		    sumBX += latQResult[i].bxValue;
		    numValid++;
		}
	    }

	    latQuality->bxValue = sumBX / numValid;
	    latQuality->quality = HIGHQ;

	    if(debug) std::cout<<"DTp2:evaluateLateralQuality \t\t\t\t\t Lateralidad ACEPTADA. HIGHQ."<<std::endl;
	}
    // Sólo una traza disjunta válida => Traza de muón incompleta pero válida.
    else { 
	if (latQResult[0].latQValid or latQResult[1].latQValid or
	    latQResult[2].latQValid or latQResult[3].latQValid)
	    {
		latQuality->valid   = true;
		latQuality->quality = LOWQ;
		for (int i = 0; i < 4; i++)
		    if (latQResult[i].latQValid) {
			latQuality->bxValue = latQResult[i].bxValue;
			/*
			 * En los casos que haya una combinación de 4 hits válidos pero
			 * sólo 3 de ellos formen traza (calidad 2), esto permite detectar
			 * la layer con el hit que no encaja en la recta, y así poder
			 * invalidarlo, cambiando su valor por "-1" como si de una mezcla
			 * de 3 hits pura se tratara.
			 * Esto es útil para los filtros posteriores.
			 */
			latQuality->invalidateHitIdx = getOmittedHit( i );
			break;
		    }

		if(debug) std::cout<<"DTp2:evaluateLateralQuality \t\t\t\t\t Lateralidad ACEPTADA. LOWQ."<<std::endl;
	    }
	else {
	    if(debug) std::cout<<"DTp2:evaluateLateralQuality \t\t\t\t\t Lateralidad DESCARTADA. NOPATH."<<std::endl;
	}
    }
}

/**
 * Valida, para una combinación de capas (3), celdas y lateralidad, si los
 * valores temporales cumplen el criterio de mean-timer.
 * En vez de comparar con un 0 estricto, que es el resultado aritmético de las
 * ecuaciones usadas de base, se incluye en la clase un valor de tolerancia
 * que por defecto vale cero, pero que se puede ajustar a un valor más
 * adecuado
 *
 * En esta primera versión de la clase, el código de generación de ecuaciones
 * se incluye en esta función, lo que es ineficiente porque obliga a calcular
 * un montón de constantes, fijas para cada combinación de celdas, que
 * tendrían que evaluarse una sóla vez en el constructor de la clase.
 * Esta disposición en el constructor estaría más proxima a la realización que
 * se tiene que llevar a término en la FPGA (en tiempo de síntesis).
 * De momento se deja aquí porque así se entiende la lógica mejor, al estar
 * descrita de manera lineal en un sólo método.
 */
void MuonPathAnalyzerPerSL::validate(LATERAL_CASES sideComb[3], int layerIndex[3],MuonPath* mPath, PARTIAL_LATQ_TYPE *latq)
{
    // Valor por defecto.
    latq->bxValue   = 0;
    latq->latQValid = false;
  
    if(debug) std::cout<<"DTp2:validate \t\t\t\t\t\t\t In validate Iniciando validacion de MuonPath para capas: "
		       <<layerIndex[0]<<"/"
		       <<layerIndex[1]<<"/"
		       <<layerIndex[2]<<std::endl;
  
    if(debug) std::cout<<"DTp2:validate \t\t\t\t\t\t\t Lateralidades parciales: "
		       <<sideComb[0]<<"/"
		       <<sideComb[1]<<"/"
		       <<sideComb[2]<<std::endl;
  
    /* Primero evaluamos si, para la combinación concreta de celdas en curso, el
       número de celdas con dato válido es 3. Si no es así, sobre esa
       combinación no se puede aplicar el mean-timer y devolvemos "false" */
    int validCells = 0;
    for (int j = 0; j < 3; j++)
	if (mPath->getPrimitive(layerIndex[j])->isValidTime()) validCells++;
  
    if (validCells != 3) {
	if(debug) std::cout<<"DTp2:validate \t\t\t\t\t\t\t No hay 3 celdas validas."<<std::endl;
	return;
    }

    if(debug) std::cout<<"DTp2:validate \t\t\t\t\t\t\t Valores de TDC: "
		       <<mPath->getPrimitive(layerIndex[0])->getTDCTime()<<"/"
		       <<mPath->getPrimitive(layerIndex[1])->getTDCTime()<<"/"
		       <<mPath->getPrimitive(layerIndex[2])->getTDCTime()<<"."
		       <<std::endl;

    if(debug) std::cout<<"DTp2:validate \t\t\t\t\t\t\t Valid TIMES: "
		       <<mPath->getPrimitive(layerIndex[0])->isValidTime()<<"/"
		       <<mPath->getPrimitive(layerIndex[1])->isValidTime()<<"/"
		       <<mPath->getPrimitive(layerIndex[2])->isValidTime()<<"."
		       <<std::endl;

  
    /* Distancias verticales entre capas inferior/media y media/superior */
    int dVertMI = layerIndex[1] - layerIndex[0];
    int dVertSM = layerIndex[2] - layerIndex[1];

    /* Distancias horizontales entre capas inferior/media y media/superior */
    int dHorzMI = cellLayout[layerIndex[1]] - cellLayout[layerIndex[0]];
    int dHorzSM = cellLayout[layerIndex[2]] - cellLayout[layerIndex[1]];

    /* Índices de pares de capas sobre las que se está actuando
       SM => Superior + Intermedia
       MI => Intermedia + Inferior
       Jugamos con los punteros para simplificar el código */
    int *layPairSM = &layerIndex[1];
    int *layPairMI = &layerIndex[0];

    /* Pares de combinaciones de celdas para composición de ecuación. Sigue la
       misma nomenclatura que el caso anterior */
    LATERAL_CASES smSides[2], miSides[2];

    /* Teniendo en cuenta que en el índice 0 de "sideComb" se almacena la
       lateralidad de la celda inferior, jugando con aritmética de punteros
       extraemos las combinaciones de lateralidad para los pares SM y MI */

    memcpy(smSides, &sideComb[1], 2 * sizeof(LATERAL_CASES));
  
    memcpy(miSides, &sideComb[0], 2 * sizeof(LATERAL_CASES));
  
    float bxValue = 0;
    int coefsAB[2] = {0, 0}, coefsCD[2] = {0, 0};
    /* It's neccesary to be careful with that pointer's indirection. We need to
       retrieve the lateral coeficientes (+-1) from the lower/middle and
       middle/upper cell's lateral combinations. They are needed to evaluate the
       existance of a possible BX value, following it's calculation equation */
    getLateralCoeficients(miSides, coefsAB);
    getLateralCoeficients(smSides, coefsCD);

    /* Cada para de sumas de los 'coefsCD' y 'coefsAB' dan siempre como resultado
       0, +-2.

       A su vez, y pese a que las ecuaciones se han construido de forma genérica
       para cualquier combinación de celdas de la cámara, los valores de 'dVertMI' y
       'dVertSM' toman valores 1 o 2 puesto que los pares de celdas con los que se
       opera en realidad, o bien están contiguos, o bien sólo están separadas por
       una fila de celdas intermedia. Esto es debido a cómo se han combinado los
       grupos de celdas, para aplicar el mean-timer, en 'LAYER_ARRANGEMENTS'.

       El resultado final es que 'denominator' es siempre un valor o nulo, o
       múltiplo de 2 */
    int denominator = dVertMI * (coefsCD[1] + coefsCD[0]) -
	dVertSM * (coefsAB[1] + coefsAB[0]);

    if(denominator == 0) {
	if(debug) std::cout<<"DTp2:validate \t\t\t\t\t\t\t Imposible calcular BX. Denominador para BX = 0."<<std::endl;
	return;
    }

    /* Esta ecuación ha de ser optimizada, especialmente en su implementación
       en FPGA. El 'denominator' toma siempre valores múltiplo de 2 o nulo, por lo
       habría que evitar el cociente y reemplazarlo por desplazamientos de bits */
    bxValue = (
	       dVertMI*(dHorzSM*MAXDRIFT + eqMainBXTerm(smSides, layPairSM, mPath)) -
	       dVertSM*(dHorzMI*MAXDRIFT + eqMainBXTerm(miSides, layPairMI, mPath))
	       ) / denominator;

    if(bxValue < 0) {
	if(debug) std::cout<<"DTp2:validate \t\t\t\t\t\t\t Combinacion no valida. BX Negativo."<<std::endl;
	return;
    }

    // Redondeo del valor del tiempo de BX al nanosegundo
    if ( (bxValue - int(bxValue)) >= 0.5 ) bxValue = float(int(bxValue + 1));
    else bxValue = float(int(bxValue));

    /* Ciertos valores del tiempo de BX, siendo positivos pero objetivamente no
       válidos, pueden dar lugar a que el discriminador de traza asociado de un
       valor aparentemente válido (menor que la tolerancia y típicamente 0). Eso es
       debido a que el valor de tiempo de BX es mayor que algunos de los tiempos
       de TDC almacenados en alguna de las respectivas 'DTPrimitives', lo que da
       lugar a que, cuando se establece el valore de BX para el 'MuonPath', se
       obtengan valores de tiempo de deriva (*NO* tiempo de TDC) en la 'DTPrimitive'
       nulos, o inconsistentes, a causa de la resta entre enteros.

       Así pues, se impone como criterio de validez adicional que el valor de tiempo
       de BX (bxValue) sea siempre superior a cualesquiera valores de tiempo de TDC
       almacenados en las 'DTPrimitives' que forman el 'MuonPath' que se está
       analizando.
       En caso contrario, se descarta como inválido */

    for (int i = 0; i < 3; i++)
	if (mPath->getPrimitive(layerIndex[i])->isValidTime()) {
	    int diffTime =
		mPath->getPrimitive(layerIndex[i])->getTDCTimeNoOffset() - bxValue;

	    if (diffTime < 0 or diffTime > MAXDRIFT) {
		if(debug) std::cout<<"DTp2:validate \t\t\t\t\t\t\t Valor de BX inválido. Al menos un tiempo de TDC sin sentido"<<std::endl;
		return;
	    }
	}

    /* Si se llega a este punto, el valor de BX y la lateralidad parcial se dan
     * por válidas.
     */
    latq->bxValue   = bxValue;
    latq->latQValid = true;
}//finish validate

/**
 * Evalúa la suma característica de cada par de celdas, según la lateralidad
 * de la trayectoria.
 * El orden de los índices de capa es crítico:
 *    layerIdx[0] -> Capa más baja,
 *    layerIdx[1] -> Capa más alta
 */
int MuonPathAnalyzerPerSL::eqMainBXTerm(LATERAL_CASES sideComb[2], int layerIdx[2],MuonPath* mPath)
{
    int eqTerm = 0, coefs[2];
    
    getLateralCoeficients(sideComb, coefs);
    
    eqTerm = coefs[0] * mPath->getPrimitive(layerIdx[0])->getTDCTimeNoOffset() +
	coefs[1] * mPath->getPrimitive(layerIdx[1])->getTDCTimeNoOffset();
    
    if(debug) std::cout<<"DTp2:eqMainBXTerm \t\t\t\t\t In eqMainBXTerm EQTerm(BX): "<<eqTerm<<std::endl;
    
    return (eqTerm);
}

/**
 * Evalúa la suma característica de cada par de celdas, según la lateralidad
 * de la trayectoria. Semejante a la anterior, pero aplica las correcciones
 * debidas a los retardos de la electrónica, junto con la del Bunch Crossing
 *
 * El orden de los índices de capa es crítico:
 *    layerIdx[0] -> Capa más baja,
 *    layerIdx[1] -> Capa más alta
 */
int MuonPathAnalyzerPerSL::eqMainTerm(LATERAL_CASES sideComb[2], int layerIdx[2],MuonPath* mPath, int bxValue)
{
    int eqTerm = 0, coefs[2];

    getLateralCoeficients(sideComb, coefs);
    
    eqTerm = coefs[0] * (mPath->getPrimitive(layerIdx[0])->getTDCTimeNoOffset() -
			 bxValue) +
	coefs[1] * (mPath->getPrimitive(layerIdx[1])->getTDCTimeNoOffset() -
		    bxValue);
    
    if(debug) std::cout<<"DTp2:\t\t\t\t\t EQTerm(Main): "<<eqTerm<<std::endl;
    
    return (eqTerm);
}

/**
 * Devuelve los coeficientes (+1 ó -1) de lateralidad para un par dado.
 * De momento es útil para poder codificar la nueva funcionalidad en la que se
 * calcula el BX.
 */

void MuonPathAnalyzerPerSL::getLateralCoeficients(LATERAL_CASES sideComb[2],int *coefs)
{
    if ((sideComb[0] == LEFT) && (sideComb[1] == LEFT)) {
	*(coefs)     = +1;
	*(coefs + 1) = -1;
    }
    else if ((sideComb[0] == LEFT) && (sideComb[1] == RIGHT)){
	*(coefs)     = +1;
	*(coefs + 1) = +1;
    }
    else if ((sideComb[0] == RIGHT) && (sideComb[1] == LEFT)){
	*(coefs)     = -1;
	*(coefs + 1) = -1;
    }
    else if ((sideComb[0] == RIGHT) && (sideComb[1] == RIGHT)){
	*(coefs)     = -1;
	*(coefs + 1) = +1;
    }
}

/**
 * Determines if all valid partial lateral combinations share the same value
 * of 'bxValue'.
 */
bool MuonPathAnalyzerPerSL::sameBXValue(PARTIAL_LATQ_TYPE* latq) {

    bool result = true;
    /*
      Para evitar los errores de precision en el cálculo, en vez de forzar un
      "igual" estricto a la hora de comparar los diferentes valores de BX, se
      obliga a que la diferencia entre pares sea menor que un cierto valor umbral.
      Para hacerlo cómodo se crean 6 booleanos que evalúan cada posible diferencia
    */
  
    if(debug) std::cout<<"Dtp2:sameBXValue bxTolerance: "<<bxTolerance<<std::endl;

    if(debug) std::cout<<"Dtp2:sameBXValue \t\t\t\t\t\t d01:"<<abs(latq[0].bxValue - latq[1].bxValue)<<std::endl;
    if(debug) std::cout<<"Dtp2:sameBXValue \t\t\t\t\t\t d02:"<<abs(latq[0].bxValue - latq[2].bxValue)<<std::endl;
    if(debug) std::cout<<"Dtp2:sameBXValue \t\t\t\t\t\t d03:"<<abs(latq[0].bxValue - latq[3].bxValue)<<std::endl;
    if(debug) std::cout<<"Dtp2:sameBXValue \t\t\t\t\t\t d12:"<<abs(latq[1].bxValue - latq[2].bxValue)<<std::endl;
    if(debug) std::cout<<"Dtp2:sameBXValue \t\t\t\t\t\t d13:"<<abs(latq[1].bxValue - latq[3].bxValue)<<std::endl;
    if(debug) std::cout<<"Dtp2:sameBXValue \t\t\t\t\t\t d23:"<<abs(latq[2].bxValue - latq[3].bxValue)<<std::endl;

    bool d01, d02, d03, d12, d13, d23;
    d01 = (abs(latq[0].bxValue - latq[1].bxValue) <= bxTolerance) ? true : false;
    d02 = (abs(latq[0].bxValue - latq[2].bxValue) <= bxTolerance) ? true : false;
    d03 = (abs(latq[0].bxValue - latq[3].bxValue) <= bxTolerance) ? true : false;
    d12 = (abs(latq[1].bxValue - latq[2].bxValue) <= bxTolerance) ? true : false;
    d13 = (abs(latq[1].bxValue - latq[3].bxValue) <= bxTolerance) ? true : false;
    d23 = (abs(latq[2].bxValue - latq[3].bxValue) <= bxTolerance) ? true : false;

    /* Casos con 4 grupos de combinaciones parciales de lateralidad validas */
    if ((latq[0].latQValid && latq[1].latQValid && latq[2].latQValid &&
	 latq[3].latQValid) && !(d01 && d12 && d23))
	result = false;
    else
	/* Los 4 casos posibles de 3 grupos de lateralidades parciales validas */
	if ( ((latq[0].latQValid && latq[1].latQValid && latq[2].latQValid) &&
	      !(d01 && d12)
	      )
	     or
	     ((latq[0].latQValid && latq[1].latQValid && latq[3].latQValid) &&
	      !(d01 && d13)
	      )
	     or
	     ((latq[0].latQValid && latq[2].latQValid && latq[3].latQValid) &&
	      !(d02 && d23)
	      )
	     or
	     ((latq[1].latQValid && latq[2].latQValid && latq[3].latQValid) &&
	      !(d12 && d23)
	      )
	     )
	    result = false;
	else
	    /* Por ultimo, los 6 casos posibles de pares de lateralidades parciales validas */

	    if ( ((latq[0].latQValid && latq[1].latQValid) && !d01) or
		 ((latq[0].latQValid && latq[2].latQValid) && !d02) or
		 ((latq[0].latQValid && latq[3].latQValid) && !d03) or
		 ((latq[1].latQValid && latq[2].latQValid) && !d12) or
		 ((latq[1].latQValid && latq[3].latQValid) && !d13) or
		 ((latq[2].latQValid && latq[3].latQValid) && !d23) )
		result = false;
  
    return result;
}

/** Calcula los parámetros de la(s) trayectoria(s) detectadas.
 *
 * Asume que el origen de coordenadas está en al lado 'izquierdo' de la cámara
 * con el eje 'X' en la posición media vertical de todas las celdas.
 * El eje 'Y' se apoya sobre los hilos de las capas 1 y 3 y sobre los costados
 * de las capas 0 y 2.
 */
void MuonPathAnalyzerPerSL::calculatePathParameters(MuonPath* mPath) {
    // El orden es importante. No cambiar sin revisar el codigo.
    if(debug) std::cout<<"DTp2:calculatePathParameters \t\t\t\t\t\t  calculating calcCellDriftAndXcoor(mPath) "<<std::endl;
    calcCellDriftAndXcoor(mPath);
    //calcTanPhiXPosChamber(mPath);
    if(debug) std::cout<<"DTp2:calculatePathParameters \t\t\t\t\t\t  checking mPath->getQuality() "<<mPath->getQuality()<<std::endl;
    if (mPath->getQuality() == HIGHQ or mPath->getQuality() == HIGHQGHOST){
	if(debug) std::cout<<"DTp2:calculatePathParameters \t\t\t\t\t\t\t  Quality test passed, now calcTanPhiXPosChamber4Hits(mPath) "<<std::endl;
	calcTanPhiXPosChamber4Hits(mPath);
    }else{
	if(debug) std::cout<<"DTp2:calculatePathParameters \t\t\t\t\t\t\t  Quality test NOT passed calcTanPhiXPosChamber3Hits(mPath) "<<std::endl;
	calcTanPhiXPosChamber3Hits(mPath);
    }
    if(debug) std::cout<<"DTp2:calculatePathParameters \t\t\t\t\t\t calcChiSquare(mPath) "<<std::endl;
    calcChiSquare(mPath);
}

void MuonPathAnalyzerPerSL::calcTanPhiXPosChamber(MuonPath* mPath)
{
    /*
      La mayoría del código de este método tiene que ser optimizado puesto que
      se hacen llamadas y cálculos redundantes que ya se han evaluado en otros
      métodos previos.

      Hay que hacer una revisión de las ecuaciones para almacenar en el 'MuonPath'
      una serie de parámetro característicos (basados en sumas y productos, para
      que su implementación en FPGA sea sencilla) con los que, al final del
      proceso, se puedan calcular el ángulo y la coordenada horizontal.

      De momento se deja este código funcional extraído directamente de las
      ecuaciones de la recta.
    */
    int layerIdx[2];
    /*
      To calculate path's angle are only necessary two valid primitives.
      This method should be called only when a 'MuonPath' is determined as valid,
      so, at least, three of its primitives must have a valid time.
      With this two comparitions (which can be implemented easily as multiplexors
      in the FPGA) this method ensures to catch two of those valid primitives to
      evaluate the angle.

      The first one is below the middle line of the superlayer, while the other
      one is above this line
    */
    if (mPath->getPrimitive(0)->isValidTime()) layerIdx[0] = 0;
    else layerIdx[0] = 1;

    if (mPath->getPrimitive(3)->isValidTime()) layerIdx[1] = 3;
    else layerIdx[1] = 2;

    /* We identify along which cells' sides the muon travels */
    LATERAL_CASES sideComb[2];
    sideComb[0] = (mPath->getLateralComb())[ layerIdx[0] ];
    sideComb[1] = (mPath->getLateralComb())[ layerIdx[1] ];

    /* Horizontal gap between cells in cell's semi-length units */
    int dHoriz = (mPath->getCellHorizontalLayout())[ layerIdx[1] ] -
	(mPath->getCellHorizontalLayout())[ layerIdx[0] ];

    /* Vertical gap between cells in cell's height units */
    int dVert = layerIdx[1] -layerIdx[0];

    /*-----------------------------------------------------------------*/
    /*--------------------- Phi angle calculation ---------------------*/
    /*-----------------------------------------------------------------*/
    float num = CELL_SEMILENGTH * dHoriz +
	DRIFT_SPEED *
	eqMainTerm(sideComb, layerIdx, mPath,
		   mPath->getBxTimeValue()
		   );

    float denom = CELL_HEIGHT * dVert;
    float tanPhi = num / denom;

    mPath->setTanPhi(tanPhi);

    /*-----------------------------------------------------------------*/
    /*----------------- Horizontal coord. calculation -----------------*/
    /*-----------------------------------------------------------------*/

    /*
      Using known coordinates, relative to superlayer axis reference, (left most
      superlayer side, and middle line between 2nd and 3rd layers), calculating
      horizontal coordinate implies using a basic line equation:
      (y - y0) = (x - x0) * cotg(Phi)
      This horizontal coordinate can be obtained setting y = 0 on last equation,
      and also setting y0 and x0 with the values of a known muon's path cell
      position hit.
      It's enough to use the lower cell (layerIdx[0]) coordinates. So:
      xC = x0 - y0 * tan(Phi)
    */
    float lowerXPHorizPos = mPath->getXCoorCell( layerIdx[0] );

    float lowerXPVertPos = 0; // This is only the absolute value distance.
    if (layerIdx[0] == 0) lowerXPVertPos = CELL_HEIGHT + CELL_SEMIHEIGHT;
    else                  lowerXPVertPos = CELL_SEMIHEIGHT;

    mPath->setHorizPos( lowerXPHorizPos + lowerXPVertPos * tanPhi );
}

/**
 * Cálculos de coordenada y ángulo para un caso de 4 HITS de alta calidad.
 */
void MuonPathAnalyzerPerSL::calcTanPhiXPosChamber4Hits(MuonPath* mPath) {
    float tanPhi = (3 * mPath->getXCoorCell(3) +
		    mPath->getXCoorCell(2) -
		    mPath->getXCoorCell(1) -
		    3 * mPath->getXCoorCell(0)) / (10 * CELL_HEIGHT);

    mPath->setTanPhi(tanPhi);

    float XPos = (mPath->getXCoorCell(0) +
		  mPath->getXCoorCell(1) +
		  mPath->getXCoorCell(2) +
		  mPath->getXCoorCell(3)) / 4;

    mPath->setHorizPos( XPos );
}

/**
 * Cálculos de coordenada y ángulo para un caso de 3 HITS.
 */
void MuonPathAnalyzerPerSL::calcTanPhiXPosChamber3Hits(MuonPath* mPath) {
    int layerIdx[2];

    if (mPath->getPrimitive(0)->isValidTime()) layerIdx[0] = 0;
    else layerIdx[0] = 1;

    if (mPath->getPrimitive(3)->isValidTime()) layerIdx[1] = 3;
    else layerIdx[1] = 2;

    /* We identify along which cells' sides the muon travels */
    LATERAL_CASES sideComb[2];
    sideComb[0] = (mPath->getLateralComb())[ layerIdx[0] ];
    sideComb[1] = (mPath->getLateralComb())[ layerIdx[1] ];

    /* Horizontal gap between cells in cell's semi-length units */
    int dHoriz = (mPath->getCellHorizontalLayout())[ layerIdx[1] ] -
	(mPath->getCellHorizontalLayout())[ layerIdx[0] ];

    /* Vertical gap between cells in cell's height units */
    int dVert = layerIdx[1] -layerIdx[0];

    /*-----------------------------------------------------------------*/
    /*--------------------- Phi angle calculation ---------------------*/
    /*-----------------------------------------------------------------*/
    float num = CELL_SEMILENGTH * dHoriz +
	DRIFT_SPEED *
	eqMainTerm(sideComb, layerIdx, mPath,
		   mPath->getBxTimeValue()
		   );

    float denom = CELL_HEIGHT * dVert;
    float tanPhi = num / denom;

    mPath->setTanPhi(tanPhi);

    /*-----------------------------------------------------------------*/
    /*----------------- Horizontal coord. calculation -----------------*/
    /*-----------------------------------------------------------------*/
    float XPos = 0;
    if (mPath->getPrimitive(0)->isValidTime() and
	mPath->getPrimitive(3)->isValidTime())
	XPos = (mPath->getXCoorCell(0) + mPath->getXCoorCell(3)) / 2;
    else
	XPos = (mPath->getXCoorCell(1) + mPath->getXCoorCell(2)) / 2;

    mPath->setHorizPos( XPos );
}

/**
 * Calcula las distancias de deriva respecto de cada "wire" y la posición
 * horizontal del punto de interacción en cada celda respecto del sistema
 * de referencia de la cámara.
 *
 * La posición horizontal de cada hilo es calculada en el "DTPrimitive".
 */
void MuonPathAnalyzerPerSL::calcCellDriftAndXcoor(MuonPath *mPath) {
    //Distancia de deriva en la celda respecto del wire". NO INCLUYE SIGNO.
    float driftDistance;
    float wireHorizPos; // Posicion horizontal del wire.
    float hitHorizPos;  // Posicion del muon en la celda.

    for (int i = 0; i <= 3; i++)
	if (mPath->getPrimitive(i)->isValidTime()) {
	    // Drift distance.
	    driftDistance = DRIFT_SPEED *
		( mPath->getPrimitive(i)->getTDCTimeNoOffset() -
		  mPath->getBxTimeValue()
		  );

	    wireHorizPos = mPath->getPrimitive(i)->getWireHorizPos();

	    if ( (mPath->getLateralComb())[ i ] == LEFT )
		hitHorizPos = wireHorizPos - driftDistance;
	    else
		hitHorizPos = wireHorizPos + driftDistance;

	    mPath->setXCoorCell(hitHorizPos, i);
	    mPath->setDriftDistance(driftDistance, i);
	}
}

/**
 * Calcula el estimador de calidad de la trayectoria.
 */
void MuonPathAnalyzerPerSL::calcChiSquare(MuonPath *mPath) {

    float xi, zi, factor;

    float chi = 0;
    float mu  = mPath->getTanPhi();
    float b   = mPath->getHorizPos();

    const float baseWireYPos = -1.5 * CELL_HEIGHT;

    for (int i = 0; i <= 3; i++)
	if ( mPath->getPrimitive(i)->isValidTime() ) {
	    zi = baseWireYPos + CELL_HEIGHT * i;
	    xi = mPath->getXCoorCell(i);

	    factor = xi - mu*zi - b;
	    chi += (factor * factor);
	}
    mPath->setChiSq(chi);
}
int MuonPathAnalyzerPerSL::getOmittedHit(int idx) {
  
    int ans = -1;
  
    switch (idx) {
    case 0: ans = 3; break;
    case 1: ans = 0; break;
    case 2: ans = 2; break;
    case 3: ans = 1; break;
    }

    return ans;
}
bool MuonPathAnalyzerPerSL::hasPosRF(int wh,int sec){
    return  wh>0 || (wh==0 && sec%4>1);
}

void MuonPathAnalyzerPerSL::analyze(MuonPath *inMPath, std::vector<metaPrimitive> &metaPrimitives){
  if(debug) std::cout<<"DTp2:analyze \t\t\t\t starts"<<std::endl;
  if(debug) std::cout<<"DTp2:analyze \t\t\t\t In analyze function checking if inMPath->isAnalyzable() "<<inMPath->isAnalyzable()<<std::endl;

    // LOCATE MPATH:
  int selected_Id=0;
  if     (inMPath->getPrimitive(0)->getTDCTime()!=-1) selected_Id=inMPath->getPrimitive(0)->getCameraId();
  else if(inMPath->getPrimitive(1)->getTDCTime()!=-1) selected_Id=inMPath->getPrimitive(1)->getCameraId(); 
  else if(inMPath->getPrimitive(2)->getTDCTime()!=-1) selected_Id=inMPath->getPrimitive(2)->getCameraId(); 
  else if(inMPath->getPrimitive(3)->getTDCTime()!=-1) selected_Id=inMPath->getPrimitive(3)->getCameraId(); 
  
  DTLayerId thisLId(selected_Id);
  if(debug) std::cout<<"Building up MuonPathSLId from rawId in the Primitive"<<std::endl;
  DTSuperLayerId MuonPathSLId(thisLId.wheel(),thisLId.station(),thisLId.sector(),thisLId.superLayer());
  if(debug) std::cout<<"The MuonPathSLId is"<<MuonPathSLId<<std::endl;
  DTWireId wireId(MuonPathSLId,2,1);
  // END LOCATE MPATH (MAYBE WE NEED TO RE-WRITE THIS)
  if (chosen_sl<4 && thisLId.superLayer()!=chosen_sl) return; // avoid running when mpath not in chosen SL (for 1SL fitting)

  // Clonamos el objeto analizado.
  MuonPath *mPath = new MuonPath(inMPath);
  
  int wi[8],tdc[8];
  DTPrimitive Prim0(mPath->getPrimitive(0)); wi[0]=Prim0.getChannelId();tdc[0]=Prim0.getTDCTime();
  DTPrimitive Prim1(mPath->getPrimitive(1)); wi[1]=Prim1.getChannelId();tdc[1]=Prim1.getTDCTime();
  DTPrimitive Prim2(mPath->getPrimitive(2)); wi[2]=Prim2.getChannelId();tdc[2]=Prim2.getTDCTime();
  DTPrimitive Prim3(mPath->getPrimitive(3)); wi[3]=Prim3.getChannelId();tdc[3]=Prim3.getTDCTime();
  for(int i=4;i<8;i++){wi[i]=-1;tdc[i]=-1;}
  
  if (mPath->isAnalyzable()) {
    if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t yes it is analyzable "<<mPath->isAnalyzable()<<std::endl;
    setCellLayout( mPath->getCellHorizontalLayout() );
    evaluatePathQuality(mPath);
  }else{
    if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t no it is NOT analyzable "<<mPath->isAnalyzable()<<std::endl;
  }
  
    
  if(debug) std::cout<<"DTp2:analyze \t\t\t\t checking if it passes the min quality cut "<<mPath->getQuality()<<">"<<minQuality<<std::endl;
  if ( mPath->getQuality() >= minQuality ) {
    if(debug) std::cout<<"DTp2:analyze \t\t\t\t min quality achievedCalidad: "<<mPath->getQuality()<<std::endl;
    for (int i = 0; i <= 3; i++){
      if(debug) std::cout<<"DTp2:analyze \t\t\t\t  Capa: "<<mPath->getPrimitive(i)->getLayerId()
			 <<" Canal: "<<mPath->getPrimitive(i)->getChannelId()<<" TDCTime: "<<mPath->getPrimitive(i)->getTDCTime()<<std::endl;
    }
    if(debug) std::cout<<"DTp2:analyze \t\t\t\t Starting lateralities loop, totalNumValLateralities: "<<totalNumValLateralities<<std::endl;
    
    float best_chi2=99999.;
    float chi2_jm_tanPhi=999;
    float chi2_jm_x=-1;
    float chi2_jm_t0=-1;
    float chi2_phi=-1;
    float chi2_phiB=-1;
    float chi2_chi2=-1;
    MP_QUALITY chi2_quality=NOPATH;
    
    for (int i = 0; i < totalNumValLateralities; i++) {//here
      if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t laterality #- "<<i<<std::endl;
      if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t laterality #- "<<i<<" checking quality:"<<std::endl;
      if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t laterality #- "<<i<<" checking mPath Quality="<<mPath->getQuality()<<std::endl;   
      if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t laterality #- "<<i<<" latQuality[i].val="<<latQuality[i].valid<<std::endl;   
      if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t laterality #- "<<i<<" before if:"<<std::endl;
      
      if (latQuality[i].valid and (((mPath->getQuality()==HIGHQ or mPath->getQuality()==HIGHQGHOST) and latQuality[i].quality==HIGHQ)
				   or
				   ((mPath->getQuality() == LOWQ or mPath->getQuality()==LOWQGHOST) and latQuality[i].quality==LOWQ))){
	
	if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t laterality #- "<<i<<" inside if"<<std::endl;
	mPath->setBxTimeValue(latQuality[i].bxValue);
	if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t laterality #- "<<i<<" settingLateralCombination"<<std::endl;
	mPath->setLateralComb(lateralities[i]);
	if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t laterality #- "<<i<<" done settingLateralCombination"<<std::endl;
	
	MuonPath *mpAux = new MuonPath(mPath);
	
	int idxHitNotValid = latQuality[i].invalidateHitIdx;
	if (idxHitNotValid >= 0) {
	  delete mpAux->getPrimitive(idxHitNotValid);
	  mpAux->setPrimitive(std::move(new DTPrimitive()), idxHitNotValid);
	}
	
	if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t  calculating parameters "<<std::endl;
	calculatePathParameters(mpAux);
	/* 
	 * Si, tras calcular los parámetros, y si se trata de un segmento
	 * con 4 hits, el chi2 resultante es superior al umbral programado,
	 * lo eliminamos y no se envía al exterior.
	 * Y pasamos al siguiente elemento.
	 */
	if ((mpAux->getQuality() == HIGHQ or mpAux->getQuality() == HIGHQGHOST) && mpAux->getChiSq() > chiSquareThreshold) {//check this if!!!
	  if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t  HIGHQ or HIGHQGHOST but min chi2 or Q test not satisfied "<<std::endl;
	}				
	else{
	  if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t  inside else, returning values: "<<std::endl;
	  if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t  BX Time = "<<mpAux->getBxTimeValue()<<std::endl;
	  if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t  BX Id   = "<<mpAux->getBxNumId()<<std::endl;
	  if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t  XCoor   = "<<mpAux->getHorizPos()<<std::endl;
	  if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t  tan(Phi)= "<<mpAux->getTanPhi()<<std::endl;
	  if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t  chi2= "<<mpAux->getChiSq()<<std::endl;
	  if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t  lateralities = "
			     <<" "<<mpAux->getLateralComb()[0]
			     <<" "<<mpAux->getLateralComb()[1]
			     <<" "<<mpAux->getLateralComb()[2]
			     <<" "<<mpAux->getLateralComb()[3]
			     <<std::endl;

	  //changing to chamber frame or reference:	  
	  DTChamberId ChId(MuonPathSLId.wheel(),MuonPathSLId.station(),MuonPathSLId.sector());

	  float jm_tanPhi=-1.*mpAux->getTanPhi(); //testing with this line
	  float jm_x=(mpAux->getHorizPos()/10.)+shiftinfo[wireId.rawId()]; 
	  //0.65 is half hight of a cell needed to go to the middle of the superlayer, here we are extrapolating with the angle of the primitive!
	  jm_x=jm_x-(zinfo[wireId.rawId()]-0.65)*jm_tanPhi; 
	  float jm_t0=mpAux->getBxTimeValue();		      
	  int quality=mpAux->getQuality();
	  
	  //computing phi and phiB
	  GlobalPoint jm_x_cmssw_global = dtGeo->chamber(ChId)->toGlobal(LocalPoint(jm_x,0.,0.));//jm_x is already extrapolated to the middle of the SL
	  int thisec = MuonPathSLId.sector();
	  if(thisec==13) thisec = 4;
	  if(thisec==14) thisec = 10;
	  float phi= jm_x_cmssw_global.phi()-0.5235988*(thisec-1);
	  float psi=atan(jm_tanPhi);
	  float phiB=hasPosRF(MuonPathSLId.wheel(),MuonPathSLId.sector()) ? psi-phi : -psi-phi ;
	  float chi2= mpAux->getChiSq();
	  
	  if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t\t\t\t  pushing back metaPrimitive at x="<<jm_x<<" tanPhi:"<<jm_tanPhi<<" t0:"<<jm_t0<<std::endl;
	  
	  //keep only the values with the best chi2 among lateralities	  
	  if(mpAux->getQuality() == HIGHQ or mpAux->getQuality() == HIGHQGHOST){ 
	    if(chi2<best_chi2){
	      chi2_jm_tanPhi=jm_tanPhi;
	      chi2_jm_x=(mpAux->getHorizPos()/10.)+shiftinfo[wireId.rawId()]; 
	      chi2_jm_x=chi2_jm_x-(zinfo[wireId.rawId()]-0.65)*chi2_jm_tanPhi; //from SL to CH
	      chi2_jm_t0=mpAux->getBxTimeValue();		      
	      chi2_phi=phi;
	      chi2_phiB=phiB;
	      chi2_chi2=chi2*0.01; //in cmssw we want cm, 1 cm^2 = 100 mm^2
	      chi2_quality= mpAux->getQuality();
	      best_chi2=mpAux->getChiSq();
	    }
	  } 
	  else {//write the metaprimitive in case no HIGHQ or HIGHQGHOST
	    if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t\t\t\t  pushing back metaprimitive no HIGHQ or HIGHQGHOST"<<std::endl;
	    metaPrimitives.push_back(metaPrimitive({thisLId.rawId(),jm_t0,jm_x,jm_tanPhi,phi,phiB,chi2,quality,
		    wi[0],tdc[0],
		    wi[1],tdc[1],
		    wi[2],tdc[2],
		    wi[3],tdc[3],
		    wi[4],tdc[4],
		    wi[5],tdc[5],
		    wi[6],tdc[6],
		    wi[7],tdc[7]
		    }));
	    
	    //Loading Info into mpAux
//	    mpAux->setChiSq(chi2*0.01);
//	    mpAux->setPhi(phi);
//	    mpAux->setPhiB(phiB);
//	    mpAux->setTanPhi(jm_tanPhi);
//	    mpAux->setHorizPos(jm_x);
//	    
//	    outMPath.push_back(std::move(mpAux));
	   
	    if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t\t\t\t  done pushing back metaprimitive no HIGHQ or HIGHQGHOST"<<std::endl;
	  }				
	}
	delete mpAux;
      }
      else{
	if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t\t\t\t  latQuality[i].valid and (((mpAux->getQuality()==HIGHQ or mpAux->getQuality()==HIGHQGHOST) and latQuality[i].quality==HIGHQ) or  ((mpAux->getQuality() == LOWQ or mpAux->getQuality()==LOWQGHOST) and latQuality[i].quality==LOWQ)) not passed"<<std::endl;
      }      
    }
    
    if(chi2_jm_tanPhi!=999){//
      if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t\t\t\t  pushing back best chi2 metaPrimitive"<<std::endl;
      metaPrimitives.push_back(metaPrimitive({thisLId.rawId(),chi2_jm_t0,chi2_jm_x,chi2_jm_tanPhi,chi2_phi,chi2_phiB,chi2_chi2,chi2_quality,
	      wi[0],tdc[0],
	      wi[1],tdc[1],
	      wi[2],tdc[2],
	      wi[3],tdc[3],
	      wi[4],tdc[4],
	      wi[5],tdc[5],
	      wi[6],tdc[6],
	      wi[7],tdc[7]
	      }));

      //Loading Info into mpAux
//      mpAux->setBxTimeValue(chi2_jm_t0);
//      mpAux->setChiSq(chi2_chi2);
//      mpAux->setPhi(chi2_phi);
//      mpAux->setPhiB(chi2_phiB);
//      mpAux->setTanPhi(chi2_jm_tanPhi);
//      mpAux->setHorizPos(chi2_jm_x);
//      mpAux->setQuality(chi2_quality);
//      
//      outMPath.push_back(std::move(mpAux));
    }
  }
  delete mPath;
  if(debug) std::cout<<"DTp2:analyze \t\t\t\t finishes"<<std::endl;
}

void MuonPathAnalyzerPerSL::analyze(MuonPath *inMPath, std::vector<MuonPath*> &outMPath){
  if(debug) std::cout<<"DTp2:analyze \t\t\t\t starts"<<std::endl;
  if(debug) std::cout<<"DTp2:analyze \t\t\t\t In analyze function checking if inMPath->isAnalyzable() "<<inMPath->isAnalyzable()<<std::endl;

    // LOCATE MPATH:
  int selected_Id=0;
  if     (inMPath->getPrimitive(0)->getTDCTime()!=-1) selected_Id=inMPath->getPrimitive(0)->getCameraId();
  else if(inMPath->getPrimitive(1)->getTDCTime()!=-1) selected_Id=inMPath->getPrimitive(1)->getCameraId(); 
  else if(inMPath->getPrimitive(2)->getTDCTime()!=-1) selected_Id=inMPath->getPrimitive(2)->getCameraId(); 
  else if(inMPath->getPrimitive(3)->getTDCTime()!=-1) selected_Id=inMPath->getPrimitive(3)->getCameraId(); 
  
  DTLayerId thisLId(selected_Id);
  if(debug) std::cout<<"Building up MuonPathSLId from rawId in the Primitive"<<std::endl;
  DTSuperLayerId MuonPathSLId(thisLId.wheel(),thisLId.station(),thisLId.sector(),thisLId.superLayer());
  if(debug) std::cout<<"The MuonPathSLId is"<<MuonPathSLId<<std::endl;
  DTWireId wireId(MuonPathSLId,2,1);
  // END LOCATE MPATH (MAYBE WE NEED TO RE-WRITE THIS)
  if (chosen_sl<4 && thisLId.superLayer()!=chosen_sl) return; // avoid running when mpath not in chosen SL (for 1SL fitting)

  // Clonamos el objeto analizado.
  MuonPath *mPath = new MuonPath(inMPath); 
 
  if (mPath->isAnalyzable()) {
    if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t yes it is analyzable "<<mPath->isAnalyzable()<<std::endl;
    setCellLayout( mPath->getCellHorizontalLayout() );
    evaluatePathQuality(mPath);
  }else{
    if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t no it is NOT analyzable "<<mPath->isAnalyzable()<<std::endl;
  }
  
    
  if(debug) std::cout<<"DTp2:analyze \t\t\t\t checking if it passes the min quality cut "<<mPath->getQuality()<<">"<<minQuality<<std::endl;
  if ( mPath->getQuality() >= minQuality ) {
    if(debug) std::cout<<"DTp2:analyze \t\t\t\t min quality achievedCalidad: "<<mPath->getQuality()<<std::endl;
    for (int i = 0; i <= 3; i++){
      if(debug) std::cout<<"DTp2:analyze \t\t\t\t  Capa: "<<mPath->getPrimitive(i)->getLayerId()
			 <<" Canal: "<<mPath->getPrimitive(i)->getChannelId()<<" TDCTime: "<<mPath->getPrimitive(i)->getTDCTime()<<std::endl;
    }
    if(debug) std::cout<<"DTp2:analyze \t\t\t\t Starting lateralities loop, totalNumValLateralities: "<<totalNumValLateralities<<std::endl;
    
    float best_chi2=99999.;
    float chi2_jm_tanPhi=999;
    float chi2_jm_x=-1;
    float chi2_jm_t0=-1;
    float chi2_phi=-1;
    float chi2_phiB=-1;
    float chi2_chi2=-1;
    MP_QUALITY chi2_quality=NOPATH;
    
    for (int i = 0; i < totalNumValLateralities; i++) {//here
      if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t laterality #- "<<i<<std::endl;
      if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t laterality #- "<<i<<" checking quality:"<<std::endl;
      if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t laterality #- "<<i<<" checking mPath Quality="<<mPath->getQuality()<<std::endl;   
      if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t laterality #- "<<i<<" latQuality[i].val="<<latQuality[i].valid<<std::endl;   
      if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t laterality #- "<<i<<" before if:"<<std::endl;
      
      if (latQuality[i].valid and (((mPath->getQuality()==HIGHQ or mPath->getQuality()==HIGHQGHOST) and latQuality[i].quality==HIGHQ)
				   or
				   ((mPath->getQuality() == LOWQ or mPath->getQuality()==LOWQGHOST) and latQuality[i].quality==LOWQ))){
	
	if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t laterality #- "<<i<<" inside if"<<std::endl;
	mPath->setBxTimeValue(latQuality[i].bxValue);
	if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t laterality #- "<<i<<" settingLateralCombination"<<std::endl;
	mPath->setLateralComb(lateralities[i]);
	if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t laterality #- "<<i<<" done settingLateralCombination"<<std::endl;
	MuonPath *mpAux = new MuonPath(mPath);

	int idxHitNotValid = latQuality[i].invalidateHitIdx;
	if (idxHitNotValid >= 0) {
	  delete mpAux->getPrimitive(idxHitNotValid);
	  mpAux->setPrimitive(std::move(new DTPrimitive()), idxHitNotValid);
	}
	
	if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t  calculating parameters "<<std::endl;
	calculatePathParameters(mpAux);
	
	/* 
	 * Si, tras calcular los parámetros, y si se trata de un segmento
	 * con 4 hits, el chi2 resultante es superior al umbral programado,
	 * lo eliminamos y no se envía al exterior.
	 * Y pasamos al siguiente elemento.
	 */
	if ((mpAux->getQuality() == HIGHQ or mpAux->getQuality() == HIGHQGHOST) && mpAux->getChiSq() > chiSquareThreshold) {//check this if!!!
	  if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t  HIGHQ or HIGHQGHOST but min chi2 or Q test not satisfied "<<std::endl;
	}				
	else{
	  if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t  inside else, returning values: "<<std::endl;
	  if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t  BX Time = "<<mpAux->getBxTimeValue()<<std::endl;
	  if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t  BX Id   = "<<mpAux->getBxNumId()<<std::endl;
	  if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t  XCoor   = "<<mpAux->getHorizPos()<<std::endl;
	  if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t  tan(Phi)= "<<mpAux->getTanPhi()<<std::endl;
	  if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t  chi2= "<<mpAux->getChiSq()<<std::endl;
	  if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t  lateralities = "
			     <<" "<<mpAux->getLateralComb()[0]
			     <<" "<<mpAux->getLateralComb()[1]
			     <<" "<<mpAux->getLateralComb()[2]
			     <<" "<<mpAux->getLateralComb()[3]
			     <<std::endl;
	  
	  //changing to chamber frame or reference:	  
	  DTChamberId ChId(MuonPathSLId.wheel(),MuonPathSLId.station(),MuonPathSLId.sector());
	  
	  float jm_tanPhi=-1.*mpAux->getTanPhi(); //testing with this line
	  float jm_x=(mpAux->getHorizPos()/10.)+shiftinfo[wireId.rawId()]; 
	  //0.65 is half hight of a cell needed to go to the middle of the superlayer, here we are extrapolating with the angle of the primitive!
	  jm_x=jm_x-(zinfo[wireId.rawId()]-0.65)*jm_tanPhi; 
	  float jm_t0=mpAux->getBxTimeValue();		      
	  MP_QUALITY quality=mpAux->getQuality(); 
	  
	  //computing phi and phiB
	  GlobalPoint jm_x_cmssw_global = dtGeo->chamber(ChId)->toGlobal(LocalPoint(jm_x,0.,0.));//jm_x is already extrapolated to the middle of the SL
	  int thisec = MuonPathSLId.sector();
	  if(thisec==13) thisec = 4;
	  if(thisec==14) thisec = 10;
	  float phi= jm_x_cmssw_global.phi()-0.5235988*(thisec-1);
	  float psi=atan(jm_tanPhi);
	  float phiB=hasPosRF(MuonPathSLId.wheel(),MuonPathSLId.sector()) ? psi-phi : -psi-phi ;
	  float chi2= mpAux->getChiSq();
	  
	  if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t\t\t\t  pushing back metaPrimitive at x="<<jm_x<<" tanPhi:"<<jm_tanPhi<<" t0:"<<jm_t0<<std::endl;
	  
	  //keep only the values with the best chi2 among lateralities	  
	  if(quality == HIGHQ or quality == HIGHQGHOST){ 
	    if(chi2<best_chi2){
	      chi2_jm_tanPhi=jm_tanPhi;
	      chi2_jm_x=(mpAux->getHorizPos()/10.)+shiftinfo[wireId.rawId()]; 
	      chi2_jm_x=chi2_jm_x-(zinfo[wireId.rawId()]-0.65)*chi2_jm_tanPhi; //from SL to CH
	      chi2_jm_t0=mpAux->getBxTimeValue();		      
	      chi2_phi=phi;
	      chi2_phiB=phiB;
	      chi2_chi2=chi2*0.01; //in cmssw we want cm, 1 cm^2 = 100 mm^2
	      chi2_quality= mpAux->getQuality();
	      best_chi2=mpAux->getChiSq();
	    }
	  } 
	  else {//write the metaprimitive in case no HIGHQ or HIGHQGHOST
	    if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t\t\t\t  pushing back metaprimitive no HIGHQ or HIGHQGHOST"<<std::endl;
	    
	    //Loading Info into mpAux
	    mpAux->setBxTimeValue(jm_t0);
 	    mpAux->setHorizPos(jm_x);
 	    mpAux->setTanPhi(jm_tanPhi);	    
 	    mpAux->setPhi(phi);
 	    mpAux->setPhiB(phiB);
 	    mpAux->setChiSq(chi2*0.01);
	    mpAux->setQuality(quality);
	    
	    outMPath.push_back(mpAux);
	   
	    if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t\t\t\t  done pushing back metaprimitive no HIGHQ or HIGHQGHOST"<<std::endl;
	  }				
	}
      }
      else{
	if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t\t\t\t  latQuality[i].valid and (((mpAux->getQuality()==HIGHQ or mpAux->getQuality()==HIGHQGHOST) and latQuality[i].quality==HIGHQ) or  ((mpAux->getQuality() == LOWQ or mpAux->getQuality()==LOWQGHOST) and latQuality[i].quality==LOWQ)) not passed"<<std::endl;
      }
    }
    
    if(chi2_jm_tanPhi!=999){//
      if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t\t\t\t  pushing back best chi2 metaPrimitive"<<std::endl;
      
      //Loading Info into mpAux
      mPath->setBxTimeValue(chi2_jm_t0);
      mPath->setChiSq(chi2_chi2);
      mPath->setPhi(chi2_phi);
      mPath->setPhiB(chi2_phiB);
      mPath->setTanPhi(chi2_jm_tanPhi);
      mPath->setHorizPos(chi2_jm_x);
      mPath->setQuality(chi2_quality);
      
      outMPath.push_back(mPath);
    }
  }
  //  delete mpAux;
  if(debug) std::cout<<"DTp2:analyze \t\t\t\t finishes"<<std::endl;
}

