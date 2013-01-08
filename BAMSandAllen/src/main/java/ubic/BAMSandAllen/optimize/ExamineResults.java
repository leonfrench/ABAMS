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
package ubic.BAMSandAllen.optimize;

import java.io.File;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes;
import ubic.BAMSandAllen.ClassSelectors.BrainRegionClassSelector;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenExpressionMatrixPair;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenNomenclaturePair;
import ubic.basecode.util.FileTools;

public class ExamineResults {
    private static Log log = LogFactory.getLog( ExamineResults.class.getName() );

    public static void main( String[] args ) throws Exception {
        // TODO Auto-generated method stub
        List<String> first = FileTools.getStringListFromFile( new File(
                "/grp/java/workspace/BAMSandAllen/data/GADataRows.1257477001069.txt" ) );
        log.info( first.size() );
        log.info( first );

        boolean squareMatrix = false;

        ConnectivityAndAllenExpressionMatrixPair expressionPair = new ConnectivityAndAllenExpressionMatrixPair(
                new BrainRegionClassSelector(), true, false, squareMatrix, Double.NaN, "NewEnergies",
                AnalyzeBAMSandAllenGenes.Direction.OUTGOING );
        // forR.writeRMatrices();

        expressionPair.slimMatrices();
        expressionPair.printDimensions();
        expressionPair.removeConnectionZeroes();

        expressionPair.printDimensions();

        log.info( "________________________________________________________" );
        log.info( "________________________________________________________" );
        Set<String> colNames = new HashSet<String>( expressionPair.getAllenDataColNames() );

        ConnectivityAndAllenNomenclaturePair x = new ConnectivityAndAllenNomenclaturePair(
                new BrainRegionClassSelector(), squareMatrix, false, colNames,
                AnalyzeBAMSandAllenGenes.Direction.OUTGOING );

        // x.sameSpace();
        x.run();
        log.info( x.getCorrelation() );

        List<String> remove = new LinkedList<String>( x.getMatrixBDataRows() );
        remove.removeAll( first );

        log.info( x.correlationReducedDataMatrix( remove ) );

        log.info( "Removing " + remove.size() + " data rows" );

        x.setMatrixBDataRows( first );

        x.printDimensions();

        log.info( x.getCorrelation() );

        x.writeImages();
        x.writeRMatrices();
        // x.test( 1000 );

    }

    /**
     * I let the GA reduce the number of genes, in this file I'm merging different runs to see if they agree, they don't
     */
    public static void main2( String[] args ) throws Exception {
        // TODO Auto-generated method stub
        List<String> first = FileTools.getStringListFromFile( new File(
                "/grp/java/workspace/BAMSandAllen/results/GA/GADataRows.279.genes.txt" ) );
        log.info( first.size() );
        log.info( first );
        List<String> second = FileTools.getStringListFromFile( new File(
                "/grp/java/workspace/BAMSandAllen/results/GA/GADataRows.144.genes.txt" ) );
        log.info( second.size() );
        log.info( second );

        second.removeAll( first );
        log.info( second.size() );
    }

}
