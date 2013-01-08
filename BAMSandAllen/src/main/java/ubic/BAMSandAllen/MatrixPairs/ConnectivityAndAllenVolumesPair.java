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

import cern.colt.list.DoubleArrayList;

import ubic.BAMSandAllen.ABAMSDataMatrix;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes;
import ubic.BAMSandAllen.SetupParameters;
import ubic.BAMSandAllen.Util;
import ubic.BAMSandAllen.AllenDataLoaders.AllenAtlasAnnotationLoader;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes.Direction;
import ubic.BAMSandAllen.ClassSelectors.BrainRegionClassSelector;
import ubic.BAMSandAllen.adjacency.VoxelVolumeAdjacency;
import ubic.BAMSandAllen.adjacency.BoxDiffAdjacency;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.math.CorrelationStats;
import ubic.basecode.math.Distance;

public class ConnectivityAndAllenVolumesPair extends ConnectivityAndAllenDataPair {
    private static Log log = LogFactory.getLog( ConnectivityAndAllenVolumesPair.class.getName() );

    public ConnectivityAndAllenVolumesPair( BrainRegionClassSelector selector, boolean squareConnectivity,
            Set<String> colNames, Direction direction ) throws Exception {
        super( selector, squareConnectivity, direction );
        try {
            AllenAtlasAnnotationLoader spaceLoader = new AllenAtlasAnnotationLoader();
            matrixB = new ABAMSDataMatrix( spaceLoader.getVolumeMatrix(), "Volume", new VoxelVolumeAdjacency() );
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
        DoubleMatrix<String, String> volumes = matrixB;
        // log.info( Arrays.asList( connectionDegrees.getRow( 0 ) ).toString() );
        // log.info( Arrays.asList( volumes.getRow( 0 ) ).toString() );
        return CorrelationStats.correl( connectionDegrees.getRow( 0 ), volumes.getRow( 0 ) );
    }

    public double getRankDegreeCorrelation() {
        DoubleMatrix<String, String> connectionDegrees = Util.columnSums( matrixA );
        DoubleMatrix<String, String> volumes = matrixB;
        log.info( "conDeg: c(" + Arrays.toString( connectionDegrees.getRow( 0 ) ) + ")" );
        log.info( "c(" + Arrays.toString( volumes.getRow( 0 ) ) + ")" );
        return Util.spearmanCorrel( connectionDegrees.getRow( 0 ), volumes.getRow( 0 ) );
    }

    /**
     * @param args
     */
    public static void main( String[] args ) throws Exception {
        // get volume data
        // AllenAtlasAnnotationLoader spaceLoader = new AllenAtlasAnnotationLoader();
        // Util.writeRTable( SetupParameters.getDataFolder() + "volume.vector.txt", spaceLoader.getVolumeMatrix() );
        // System.exit( 1 );

        log.info( "Setting cols to same of the Expression Matrix" );
        boolean squareMatrix = false;
        Set<String> colNames = null;
        Direction direction = AnalyzeBAMSandAllenGenes.Direction.APPENDED;

        colNames = new HashSet<String>( ExpressionMatrixPairFactory.getUsedCols( squareMatrix, direction ) );

        ConnectivityAndAllenVolumesPair x = new ConnectivityAndAllenVolumesPair( new BrainRegionClassSelector(),
                squareMatrix, colNames, direction );
        x.run();
        // x.runAllenStyle();
        log.info( "getDegreeCorrelation:" + x.getDegreeCorrelation() );
        log.info( "Pval:" + CorrelationStats.pvalue( x.getDegreeCorrelation(), x.getMatrixA().columns() ) );

        log.info( "getRankDegreeCorrelation:" + x.getRankDegreeCorrelation() );
        log.info( "Pval:" + CorrelationStats.spearmanPvalue( x.getRankDegreeCorrelation(), x.getMatrixA().columns() ) );

    }
}
