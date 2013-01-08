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

import java.util.LinkedList;
import java.util.List;
import java.util.Random;

import javax.media.j3d.Link;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes;
import ubic.BAMSandAllen.RankedGeneListLoader;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes.Direction;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenExpressionMatrixPair;
import ubic.BAMSandAllen.MatrixPairs.ExpressionMatrixPairFactory;
import ubic.BAMSandAllen.MatrixPairs.MatrixPair;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenPartialExpressionMatrixPair.RegressMatrix;

public class TestConnectivityRobustness {
    private static Log log = LogFactory.getLog( MatrixPair.class.getName() );

    /**
     * @param args
     */
    public static void main( String[] args ) throws Exception {

        main2();
        System.exit( 1 );

        Direction direction = AnalyzeBAMSandAllenGenes.Direction.OUTGOING;
        boolean makeVirtualRegions = true;

        boolean squareMatrix = false;
        boolean logDistance = true;
        boolean slow = false;
        RegressMatrix regressType = RegressMatrix.BOTH;
        boolean useVirtual = true;
        boolean removeNonExp = true;

        ConnectivityAndAllenExpressionMatrixPair pair = ExpressionMatrixPairFactory.connectivityPartial( direction,
                slow, regressType, useVirtual, removeNonExp, logDistance );
        pair.getCorrelation();
        // ConnectivityAndAllenExpressionMatrixPair pair = ExpressionMatrixPairFactory.connectivityAndExpression(
        // direction, useVirtual, removeNonExp );

        String filename = "/grp/java/workspace/BAMSandAllen/data/rankedGenes/near final ammon/";
        filename += "LOOGenesInOrder.in.partialcon.ammon.txt.452.0.01985.topGenes.txt";
        // filename += "LOOGenesInOrder.out.partialcon.ammon.txt.374.0.016424.topGenes.txt";

        RankedGeneListLoader aLook = new RankedGeneListLoader( filename );
        pair.setMatrixBDataRows( aLook.getLines() );
        log.info( "Original:" + pair.getCorrelation() );

        log.info( pair.runJackKnife( true ) );

        // second try

        List<String> connectionRows = pair.getMatrixA().getRowNames();
        int seed = 7;
        Random r = new Random( seed );
        List<String> removeRows = new LinkedList<String>();
        for ( String region : connectionRows ) {
            boolean remove = r.nextBoolean();
            if ( remove ) {
                removeRows.add( region );
            }
        }
        pair.removeMatrixADataRows( removeRows );
        pair.printConnectionInfo();
        pair.printDimensions();
        log.info( "After remvoing half:" + pair.getCorrelation() + ", seed=" + seed );

    }

    public static void main2() throws Exception {
        Direction direction = AnalyzeBAMSandAllenGenes.Direction.INCOMING;
        boolean makeVirtualRegions = true;

        boolean squareMatrix = false;
        boolean logDistance = true;
        boolean slow = false;
        RegressMatrix regressType = RegressMatrix.BOTH;
        boolean useVirtual = true;
        boolean removeNonExp = true;
        boolean run = false;

        List<Double> resultCor = new LinkedList<Double>();
        List<Double> resultSig = new LinkedList<Double>();

        int seed = 1;
        int splits = 2;

        for ( int i = 0; i < splits; i++ ) {
            ConnectivityAndAllenExpressionMatrixPair pair = ExpressionMatrixPairFactory.connectivityPartial( direction,
                    slow, regressType, useVirtual, removeNonExp, logDistance, squareMatrix, run );

            pair.removeZeroConnectionRows();

//            String filename = "/grp/java/workspace/BAMSandAllen/data/rankedGenes/near final ammon/";
//             filename += "LOOGenesInOrder.in.partialcon.ammon.txt.452.0.01985.topGenes.txt";
//            //filename += "LOOGenesInOrder.out.partialcon.ammon.txt.374.0.016424.topGenes.txt";
//            RankedGeneListLoader aLook = new RankedGeneListLoader( filename );
//            pair.setMatrixBDataRows( aLook.getLines() );

            // remove rows
            List<String> connectionRows = pair.getMatrixA().getRowNames();
            List<String> removeRows = new LinkedList<String>();
            Random r = new Random( seed );
            for ( String region : connectionRows ) {
                int remove = ( int ) ( r.nextDouble() * splits );
                if ( remove != i ) {
                    removeRows.add( region );
                }
            }
            pair.removeMatrixADataRows( removeRows );
            pair.printConnectionInfo();
            pair.printDimensions();
            pair.run();
            double correlation = pair.getCorrelation();
            log.info( correlation );
            resultCor.add( correlation );
            resultSig.add( pair.test( 1000 ) );
        }
        log.info( "Correlations:" + resultCor );
        log.info( "Pvalues:" + resultSig );
    }
}
