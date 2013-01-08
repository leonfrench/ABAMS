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
import java.util.List;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.poi.hssf.record.formula.ExpPtg;

import ubic.BAMSandAllen.ABAMSDataMatrix;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes;
import ubic.BAMSandAllen.SetupParameters;
import ubic.BAMSandAllen.Util;
import ubic.BAMSandAllen.AllenDataLoaders.AllenAtlasAnnotationLoader;
import ubic.BAMSandAllen.AllenDataLoaders.AllenCatalogMatrices2;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes.Direction;
import ubic.BAMSandAllen.ClassSelectors.BrainRegionClassSelector;
import ubic.BAMSandAllen.adjacency.AdjacencyCompute;
import ubic.BAMSandAllen.adjacency.CorrelationAdjacency;
import ubic.BAMSandAllen.adjacency.EuclidAdjacency;
import ubic.BAMSandAllen.adjacency.LogEuclidAdjacency;
import ubic.BAMSandAllen.adjacency.VoxelVolumeAdjacency;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;

public class AllenMatrixPairFactory {
    private static Log log = LogFactory.getLog( AllenMatrixPairFactory.class.getName() );

    public static ConnectivityAndAllenSpacePair getConnectivityAndDistancePair( Direction direction,
            boolean virtualRegions, boolean logDistance ) throws Exception {
        boolean squareMatrix = false;

        Set<String> colNames = new HashSet<String>( ExpressionMatrixPairFactory.getUsedCols( squareMatrix, direction ) );

        ConnectivityAndAllenSpacePair x = new ConnectivityAndAllenSpacePair( new BrainRegionClassSelector(),
                squareMatrix, colNames, direction, logDistance );
        if ( virtualRegions ) x.makeVirtualRegions();
        x.run();
        return x;
    }

    public static AllenMatrixPair getVolumeAndExpressionPair( Direction direction, boolean removeNonExp )
            throws Exception {
        AllenAtlasAnnotationLoader spaceLoader = new AllenAtlasAnnotationLoader();
        DoubleMatrix<String, String> volMatrix = spaceLoader.getVolumeMatrix();
        ABAMSDataMatrix volABAMS = new ABAMSDataMatrix( volMatrix, "Volume", new VoxelVolumeAdjacency() );

        ABAMSDataMatrix energyABAMS = ExpressionMatrixPairFactory.getEnergyMatrix( direction, removeNonExp );

        boolean squareMatrix = false;
        Set<String> colNames = new HashSet<String>( ExpressionMatrixPairFactory.getUsedCols( squareMatrix, direction ) );

        log.info( "_____________________________" );

        AllenMatrixPair matrixPair;

        matrixPair = new AllenMatrixPair( volABAMS, energyABAMS, colNames );
        matrixPair.run();

        log.info( matrixPair.getCorrelation() );
        return matrixPair;
    }

    public static AllenMatrixPair getSpaceAndExpressionPair( Direction direction, boolean removeNonExp,
            boolean logDistance ) throws Exception {
        boolean squareMatrix = false;

        AllenAtlasAnnotationLoader spaceLoader = new AllenAtlasAnnotationLoader();
        DoubleMatrix<String, String> spaceMatrix = spaceLoader.getCenterMatrix();

        AdjacencyCompute compute;
        if ( logDistance ) {
            compute = new LogEuclidAdjacency();
        } else {
            compute = new EuclidAdjacency();
        }

        ABAMSDataMatrix spaceABAMS = new ABAMSDataMatrix( spaceMatrix, "Space", compute );

        ABAMSDataMatrix energyABAMS = ExpressionMatrixPairFactory.getEnergyMatrix( direction, removeNonExp );

        Set<String> colNames = new HashSet<String>( ExpressionMatrixPairFactory.getUsedCols( squareMatrix, direction ) );

        log.info( "_____________________________" );

        AllenMatrixPair matrixPair;

        // exp to space
        matrixPair = new AllenMatrixPair( spaceABAMS, energyABAMS, colNames );
        matrixPair.run();

        log.info( matrixPair.getCorrelation() );
        return matrixPair;
    }

    public static void main( String args[] ) throws Exception {
        boolean logDistance = true;
        boolean removeNonExp = true;
        boolean virtualRegions = true;

        Direction direction = AnalyzeBAMSandAllenGenes.Direction.INCOMING;

        //MatrixPair pair = getSpaceAndExpressionPair( Direction.INCOMING, removeNonExp, logDistance );
        // MatrixPair space = getVolumeAndExpressionPair( Direction.INCOMING, removeNonExp );
        ConnectivityAndAllenSpacePair pair = getConnectivityAndDistancePair( direction, virtualRegions, logDistance );
        pair.getCorrelation();
        pair.test( 1000 );

    }

}
