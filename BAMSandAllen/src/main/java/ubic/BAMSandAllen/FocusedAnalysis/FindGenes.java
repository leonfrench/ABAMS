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

import java.util.HashSet;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes;
import ubic.BAMSandAllen.RankedGeneListLoader;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes.Direction;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenPartialExpressionMatrixPair;
import ubic.BAMSandAllen.MatrixPairs.ExpressionMatrixPairFactory;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenPartialExpressionMatrixPair.RegressMatrix;

public class FindGenes {
    private static Log log = LogFactory.getLog( FindGenes.class.getName() );

    /**
     * @param args
     */
    public static void main( String[] args ) throws Exception {
        Direction direction = AnalyzeBAMSandAllenGenes.Direction.OUTGOING;
        boolean makeVirtualRegions = true;
        String printResult = "";
        boolean squareMatrix = true;
        boolean logDistance = true;
        boolean slow = false;
        RegressMatrix regressType = RegressMatrix.BOTH;
        boolean useVirtual = true;
        boolean removeNonExp = true;
        boolean run = false;
        ConnectivityAndAllenPartialExpressionMatrixPair pair;

        Set<String> divisions = new HashSet<String>();
        divisions.add( "Hindbrain" );
        divisions.add( "Interbrain" );
        divisions.add( "Midbrain" );
        divisions.add( "Cerebrum" );
        divisions.add( "Cerebellum" );

        String filename = "/grp/java/workspace/BAMSandAllen/data/rankedGenes/near final ammon/";
        // filename += "LOOGenesInOrder.in.partialcon.ammon.txt.452.0.01985.topGenes.txt";
        // filename += "LOOGenesInOrder.space.ammon.txt.420.0.018435.topGenes.txt";
        filename += "LOOGenesInOrder.out.partialcon.ammon.txt.374.0.016424.topGenes.txt";

        pair = ExpressionMatrixPairFactory.connectivityPartial( direction, slow, regressType, useVirtual, removeNonExp,
                logDistance, squareMatrix, run );
        pair.run();
        System.out.println( findGenes( filename, pair, 100 ) );
        System.exit( 1 );

        for ( String division : divisions ) {
            pair = ExpressionMatrixPairFactory.connectivityPartial( direction, slow, regressType, useVirtual,
                    removeNonExp, logDistance, squareMatrix, run );

            // used for direct measures:
            // squareMatrix = true; // temp useVirtual = false;// temp
            // ConnectivityAndAllenExpressionMatrixPair pair = ExpressionMatrixPairFactory.connectivityAndExpression(
            // direction, useVirtual, removeNonExp, run, squareMatrix );
            // pair.getMatrixA().setAdjacencyCompute( new IdentityAdjacency( pair.getMatrixA() ) );
            // pair.setDivision( division ); // set it to identity - temp
            // pair.getMatrixA().setAdjacencyCompute( new IdentityAdjacency( pair.getMatrixA() ) );

            pair.setDivision( division );
            log.info( "Division:" + pair.getCorrelation() );
            // pair.test( 1000 );

            printResult += "\ndivision:" + division + "\n";
            printResult += findGenes( filename, pair, 20 );

        }
        System.out.println( printResult );
    }

    private static String findGenes( String topGenesFile, ConnectivityAndAllenPartialExpressionMatrixPair pair,
            int howMany ) throws Exception {
        String printResult = "";

        RankedGeneListLoader aLook = new RankedGeneListLoader( topGenesFile );
        pair.setMatrixBDataRows( aLook.getLines().subList( 0, howMany ) );
        log.info( "Top twenty:" + pair.getCorrelation() );
        pair.test( 1000 );
        boolean add = true;

        for ( String rowName : pair.getMatrixBDataRows() ) {
            log.info( rowName + ":" );
            pair.setToSingleGene( rowName, add );
            double cor = pair.getCorrelation();
            double pval = pair.test( 1000, false );
            if ( pval < 0.05 && !Double.isNaN( cor ) ) {
                log.info( "LOW:" + rowName + "  correlation:" + cor );
                printResult += rowName + "  correlation:" + cor + " pval:" + pval + "\n";
            }
        }
        return printResult;
    }
}
