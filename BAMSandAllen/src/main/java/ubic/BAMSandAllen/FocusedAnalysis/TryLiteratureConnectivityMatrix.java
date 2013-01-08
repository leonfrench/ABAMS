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

import ubic.BAMSandAllen.SetupParameters;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenExpressionMatrixPair;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenPartialExpressionMatrixPair;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenSpacePair;
import ubic.BAMSandAllen.MatrixPairs.ExpressionMatrixPairFactory;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenPartialExpressionMatrixPair.RegressMatrix;
import ubic.BAMSandAllen.geneFilters.GeneFilter;
import ubic.BAMSandAllen.geneFilters.NonExpFilter;

public class TryLiteratureConnectivityMatrix {
    protected static Log log = LogFactory.getLog( TryLiteratureConnectivityMatrix.class );

    /**
     * @param args
     */
    public static void main( String[] args ) throws Exception {
        String base = SetupParameters.config.getString( "whitetext.iteractions.matricesFolder" );
        String filename = base;
        // filename += "Positives.rat.WhiteTextUnseen.matrix.txt.propigated";
        // filename += "Positives.WhiteTextUnseen.matrix.txt.propigated";
        // filename += "Positives.rat.WhiteTextUnseen.matrix.txt";
        // filename += "Negatives.rat.WhiteTextUnseen.matrix.txt";

        filename += "Positives.WhiteTextUnseen.matrix.txt.propigated"; // 115 columns used
        // filename += "Positives.rat.WhiteTextUnseen.matrix.txt.propigated"; //106 columns

        runProximity( filename );
        // runNormal(filename);
    }

    public static void runProximity( String filename ) throws Exception {
        RegressMatrix regressType = RegressMatrix.BOTH;
        ConnectivityAndAllenPartialExpressionMatrixPair pair = getPartialLiteratureExpressionPair( filename, regressType );

        log.info( "Correlation:" + pair.getCorrelation() );

        pair.test( 1000 );
    }

    public static ConnectivityAndAllenPartialExpressionMatrixPair getPartialLiteratureExpressionPair( String filename, RegressMatrix regressType )
            throws Exception {
        double zeroReplacement = Double.NaN;
        boolean onePlusLog = false;
        boolean doLog = true;
        String matrixName = "NewEnergies";
         

        ConnectivityAndAllenSpacePair spacePair = new ConnectivityAndAllenSpacePair( filename, null, doLog );
        spacePair.makeVirtualRegions();
        spacePair.run();

        ConnectivityAndAllenPartialExpressionMatrixPair pair = new ConnectivityAndAllenPartialExpressionMatrixPair(
                filename, doLog, onePlusLog, zeroReplacement, matrixName, spacePair.getMatrixB(), regressType );

        pair.makeVirtualRegions();

        pair.applyGeneFilter( new NonExpFilter() );

        pair.run();

        // default filters
        for ( GeneFilter filter : ExpressionMatrixPairFactory.getDefaultExpressionFilters() ) {
            pair.applyGeneFilter( filter );
        }
        return pair;
    }

    public static void runNormal( String filename ) throws Exception {
        double zeroReplacement = Double.NaN;
        boolean onePlusLog = false;
        boolean doLog = true;
        String matrixName = "NewEnergies";

        // WARNING, removes occurance counts!!
        ConnectivityAndAllenExpressionMatrixPair pair = new ConnectivityAndAllenExpressionMatrixPair( filename, doLog,
                onePlusLog, zeroReplacement, matrixName );

        pair.makeVirtualRegions();

        pair.applyGeneFilter( new NonExpFilter() );

        pair.run();

        // default filters
        for ( GeneFilter filter : ExpressionMatrixPairFactory.getDefaultExpressionFilters() ) {
            pair.applyGeneFilter( filter );
        }

        log.info( pair.getCorrelation() );

        pair.test( 1000, false );

        pair.writeImages();

    }
}
