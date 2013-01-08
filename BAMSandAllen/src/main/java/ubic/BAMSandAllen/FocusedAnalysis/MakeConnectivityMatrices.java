/*******************************************************************************
 * The ABAMS project
 * 
 * Copyright (c) 2012 University of British Columbia
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *       http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 ******************************************************************************/
package ubic.BAMSandAllen.FocusedAnalysis;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.poi.hssf.record.formula.ExpPtg;

import ubic.BAMSandAllen.ABAMSDataMatrix;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes;
import ubic.BAMSandAllen.SetupParameters;
import ubic.BAMSandAllen.StructureCatalogAnalyze;
import ubic.BAMSandAllen.Util;
import ubic.BAMSandAllen.AllenDataLoaders.StructureCatalogLoader;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes.Direction;
import ubic.BAMSandAllen.BAMSDataLoaders.BAMSDataLoader;
import ubic.BAMSandAllen.ClassSelectors.BrainRegionClassSelector;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenExpressionMatrixPair;
import ubic.BAMSandAllen.MatrixPairs.ExpressionMatrixPairFactory;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenPartialExpressionMatrixPair.RegressMatrix;
import ubic.BAMSandAllen.adjacency.DotProductAdjacency;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;

/**
 * Writes connectivity matrices for supplement website.
 * 
 * @author leon
 */
public class MakeConnectivityMatrices {
    private static Log log = LogFactory.getLog( SquareTests.class.getName() );

    public static void main( String[] args ) throws Exception {
        // ABA square (must be incoming)
        Direction direction = AnalyzeBAMSandAllenGenes.Direction.ANYDIRECTION;
        ConnectivityAndAllenExpressionMatrixPair pair;

        boolean useVirtual = true;
        boolean removeNonExp = true;
        boolean run = true;
        boolean squareMatrix = true;

        // 1155 connections
        direction = AnalyzeBAMSandAllenGenes.Direction.INCOMING; // incoming is the traditional matrix setup

        // pair = ExpressionMatrixPairFactory.connectivityAndExpression( direction, useVirtual, removeNonExp, run,
        // squareMatrix );
        // pair.printConnectionInfo();
        // pair.printDimensions();
        // 127 x 127

        // 141x141 Allen regions, some rows and columbs will be all zeroes - no reported connections
        // 1157 connections
        // pair = SquareTests.make142SquarePair( direction );
        // pair.printConnectionInfo();
        // pair.printDimensions();
        // pair.writeRMatrices();

        // ABA rectangle incoming
        // ABA recangle outgoing
        // direction = AnalyzeBAMSandAllenGenes.Direction.OUTGOING;
        // pair = ExpressionMatrixPairFactory.connectivityAndExpression( direction, useVirtual, removeNonExp );
        // pair.printConnectionInfo();
        // pair.printDimensions();
//        pair.writeRMatrices();

        direction = AnalyzeBAMSandAllenGenes.Direction.INCOMING;

        // BAMS 962 x 962
        boolean propigated = false;
        StructureCatalogAnalyze forMatrix = new StructureCatalogAnalyze( new BrainRegionClassSelector() );
        String filename = "Propigated.rdf";
        if ( !propigated ) filename = "Non" + filename;

        forMatrix.readModel( SetupParameters.getDataFolder() + filename );
        DoubleMatrix<String, String> dataMatrix = forMatrix.makeConnectionMatrix( direction );
        Util.writeRTable( SetupParameters.getDataFolder() + filename + ".matrix.txt", dataMatrix );

        // Heirarchy matrices
        StructureCatalogLoader loader = new StructureCatalogLoader();
        ABAMSDataMatrix matrix = new ABAMSDataMatrix( loader.getNomenclatureMatrix(), "Nomenclature", new DotProductAdjacency() );
//        log.info( matrix.getDimensionString());
//        Util.writeRTable( SetupParameters.getDataFolder() + "Nomenclature.matrix.Allen.txt", matrix );
//        Util.writeImage( SetupParameters.getDataFolder() + "Nomenclature.matrix.Allen.png", matrix );

        BAMSDataLoader bamsLoader = new BAMSDataLoader();
        matrix = new ABAMSDataMatrix( bamsLoader.getNomenclatureMatrix(), "Nomenclature", new DotProductAdjacency() );
        log.info( matrix.getDimensionString());
        Util.writeRTable( SetupParameters.getDataFolder() + "Nomenclature.matrix.BAMS.txt", matrix );
        Util.writeImage( SetupParameters.getDataFolder() + "Nomenclature.matrix.BAMS.png", matrix );
        
        
        // propigated and not
    }

}
