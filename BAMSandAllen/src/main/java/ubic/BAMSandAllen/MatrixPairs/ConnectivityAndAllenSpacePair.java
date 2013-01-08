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

import java.util.HashSet;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.ABAMSDataMatrix;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes;
import ubic.BAMSandAllen.AllenDataLoaders.AllenAtlasAnnotationLoader;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes.Direction;
import ubic.BAMSandAllen.ClassSelectors.BrainRegionClassSelector;
import ubic.BAMSandAllen.adjacency.AdjacencyCompute;
import ubic.BAMSandAllen.adjacency.EuclidAdjacency;
import ubic.BAMSandAllen.adjacency.LogEuclidAdjacency;

public class ConnectivityAndAllenSpacePair extends ConnectivityAndAllenDataPair {
    private static Log log = LogFactory.getLog( ConnectivityAndAllenSpacePair.class.getName() );

    public ConnectivityAndAllenSpacePair( BrainRegionClassSelector selector, boolean squareConnectivity,
            Set<String> colNames, Direction direction ) throws Exception {
        this( selector, squareConnectivity, colNames, direction, false );
    }

    public ConnectivityAndAllenSpacePair( String connectivityMatrix, Set<String> colNames, boolean logDistance )
            throws Exception {
        super( connectivityMatrix );
        setSpaceMatrix( colNames, logDistance );
    }

    public ConnectivityAndAllenSpacePair( BrainRegionClassSelector selector, boolean squareConnectivity,
            Set<String> colNames, Direction direction, boolean logDistance ) throws Exception {
        super( selector, squareConnectivity, direction );
        setSpaceMatrix( colNames, logDistance );
    }

    private void setSpaceMatrix( Set<String> colNames, boolean logDistance ) {
        try {
            AllenAtlasAnnotationLoader spaceLoader = new AllenAtlasAnnotationLoader();
            String matrixBname = "EuclidianDistance";

            AdjacencyCompute adjacencyCompute = null;
            if ( logDistance ) {
                adjacencyCompute = new LogEuclidAdjacency();
                matrixBname = "Log" + matrixBname;
            } else {
                adjacencyCompute = new EuclidAdjacency();
            }
            matrixB = new ABAMSDataMatrix( spaceLoader.getCenterMatrix(), matrixBname, adjacencyCompute );

            // reduce to the cols we were given
            if ( colNames != null ) {
                matrixB = matrixB.retainColumns( colNames );
            }
        } catch ( Exception e ) {
            e.printStackTrace();
            System.exit( 1 );
        }
    }

    /**
     * @param args
     */
    public static void main( String[] args ) throws Exception {
        log.info( "Setting cols to same of the Expression Matrix" );
        boolean squareMatrix = false;
        Direction direction = AnalyzeBAMSandAllenGenes.Direction.INCOMING;

        Set<String> colNames = new HashSet<String>( ExpressionMatrixPairFactory.getUsedCols( squareMatrix, direction ) );

        ConnectivityAndAllenSpacePair x = new ConnectivityAndAllenSpacePair( new BrainRegionClassSelector(),
                squareMatrix, colNames, direction );
        // !! using virtual regions!
        x.makeVirtualRegions();
        x.run();
        x.runAllenStyle();
        // x.writeImages();
        // x.writeRMatrices();
        // x.test( 1000, false );

    }
}
