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
package ubic.BAMSandAllen.MatrixPairs;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.ABAMSDataMatrix;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes;
import ubic.BAMSandAllen.Util;
import ubic.BAMSandAllen.AllenDataLoaders.AllenAtlasAnnotationLoader;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes.Direction;
import ubic.BAMSandAllen.ClassSelectors.BrainRegionClassSelector;
import ubic.BAMSandAllen.adjacency.BoxDiffAdjacency;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.math.CorrelationStats;

public class ConnectivityAndAllenDimensionsPair extends ConnectivityAndAllenDataPair {
    private static Log log = LogFactory.getLog( ConnectivityAndAllenDimensionsPair.class.getName() );

    public ConnectivityAndAllenDimensionsPair( BrainRegionClassSelector selector, boolean squareConnectivity,
            Set<String> colNames, Direction direction ) throws Exception {
        super( selector, squareConnectivity, direction );
        try {
            AllenAtlasAnnotationLoader spaceLoader = new AllenAtlasAnnotationLoader();
            matrixB = new ABAMSDataMatrix( spaceLoader.getDimensionsMatrix(), "Dimensions", new BoxDiffAdjacency() );
            // reduce to the cols we were given
            if ( colNames != null ) {
                matrixB = matrixB.retainColumns( colNames );
            }

        } catch ( Exception e ) {
            e.printStackTrace();
            System.exit( 1 );
        }
    }

    public double getDegreeCorrelation() {
        DoubleMatrix<String, String> connectionDegrees = Util.columnSums( matrixA );
        BoxDiffAdjacency volAdj = new BoxDiffAdjacency();
        DoubleMatrix<String, String> volumes = volAdj.getVolumes( matrixB );
        log.info( Arrays.asList( connectionDegrees.getRow( 0 ) ).toString() );
        log.info( Arrays.asList( volumes.getRow( 0 ) ).toString() );
        return CorrelationStats.correl( connectionDegrees.getRow( 0 ), volumes.getRow( 0 ) );
    }

    /**
     * @param args
     */
    public static void main( String[] args ) throws Exception {
        log.info( "Setting cols to same of the Expression Matrix" );
        boolean squareMatrix = false;
        Set<String> colNames = null;

        colNames = new HashSet<String>( ExpressionMatrixPairFactory.getUsedCols( squareMatrix,
                AnalyzeBAMSandAllenGenes.Direction.OUTGOING ) );

        ConnectivityAndAllenDimensionsPair x = new ConnectivityAndAllenDimensionsPair( new BrainRegionClassSelector(),
                squareMatrix, colNames, AnalyzeBAMSandAllenGenes.Direction.OUTGOING );
        x.run();
        x.runAllenStyle();
        log.info( "getDegreeCorrelation:" + x.getDegreeCorrelation() );

        // x.removeBedNucleiStria();

        // log.info( "After:" + x.getCorrelation( false ) );
        // x.printDimensions();
        // log.info( "Diff:" + x.bedStriaDiff() );
        // log.info( "Diff:" + x.bedStriaDiff() );
        // log.info( "Diff:" + x.bedStriaDiff() );
        // System.exit( 1 );

    }

}
