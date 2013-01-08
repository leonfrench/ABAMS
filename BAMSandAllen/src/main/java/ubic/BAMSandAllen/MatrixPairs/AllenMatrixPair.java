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

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.ABAMSDataMatrix;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes;
import ubic.BAMSandAllen.AllenDataLoaders.AllenAtlasAnnotationLoader;
import ubic.BAMSandAllen.AllenDataLoaders.StructureCatalogLoader;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes.Direction;
import ubic.BAMSandAllen.FocusedAnalysis.ExploreRegionNames;
import ubic.BAMSandAllen.adjacency.BoxDiffAdjacency;
import ubic.BAMSandAllen.adjacency.CorrelationAdjacency;
import ubic.BAMSandAllen.adjacency.DotProductAdjacency;
import ubic.BAMSandAllen.adjacency.LogEuclidAdjacency;
import ubic.BAMSandAllen.adjacency.VoxelVolumeAdjacency;
import ubic.BAMSandAllen.geneFilters.GeneFilter;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;

public class AllenMatrixPair extends MatrixPair {
    private static Log log = LogFactory.getLog( AllenMatrixPair.class.getName() );

    // public AllenMatrixPair( DoubleMatrix<String, String> matrixA, DoubleMatrix<String, String> matrixB,
    // String matrixNameA, String matrixNameB ) {
    // this( matrixA, matrixB, matrixNameA, matrixNameB, null );
    // }
    public AllenMatrixPair( ABAMSDataMatrix matrixA, ABAMSDataMatrix matrixB ) {
        this( matrixA, matrixB, null );
    }

    public AllenMatrixPair( ABAMSDataMatrix matrixA, ABAMSDataMatrix matrixB, Set<String> colNames ) {
        super();
        isInSameSpace = true;
        this.matrixA = matrixA;
        this.matrixB = matrixB;
        if ( colNames != null ) {
            this.matrixA = matrixA.retainColumns( colNames );
            this.matrixB = matrixB.retainColumns( colNames );
        }
    }

    public Set<String> convertANametoB( String aName ) {
        Set<String> result = new HashSet<String>();
        result.add( aName );
        return result;
    }


    public Set<String> convertBNametoA( String bName ) {
        Set<String> result = new HashSet<String>();
        result.add( bName );
        return result;
    }

    public double run() throws Exception {
        sameSpace();
        double correl = getCorrelation( false );
        printDimensions();
        log.info( matrixA.getName() + " and " + matrixB.getName() + " correlation is " + correl );
        testBoth( 1000, false );
        // writeImages();
        return correl;
    }

    /**
     * @param args
     */
    public static void main( String[] args ) throws Exception {

        boolean logDistance = true;
        boolean removeNonExp = true;

        Direction direction = AnalyzeBAMSandAllenGenes.Direction.OUTGOING;

        DoubleMatrix<String, String> energyMatrix = ExpressionMatrixPairFactory.getEnergyMatrix( direction,
                removeNonExp );
        ConnectivityAndAllenExpressionMatrixPair forR = ExpressionMatrixPairFactory.connectivityAndExpression(
                direction, true, removeNonExp, false, false );
        ExploreRegionNames explore = new ExploreRegionNames( forR );
        Collection<String> regions = explore.loader();
        log.info( regions.size() );

        ABAMSDataMatrix energyABAMS = ExpressionMatrixPairFactory.getEnergyMatrix( direction, removeNonExp );

        // new ABAMSDataMatrix( energyMatrix, "NewEnergies", new CorrelationAdjacency(
        // energyMatrix ) );

        AllenAtlasAnnotationLoader spaceLoader = new AllenAtlasAnnotationLoader();
        DoubleMatrix<String, String> spaceMatrix = spaceLoader.getCenterMatrix();
        ABAMSDataMatrix spaceABAMS = new ABAMSDataMatrix( spaceMatrix, "Space", new LogEuclidAdjacency() );

        ABAMSDataMatrix volABAMS = new ABAMSDataMatrix( spaceLoader.getVolumeMatrix(), "Volume",
                new VoxelVolumeAdjacency() );

        DoubleMatrix<String, String> dimsMatrix = spaceLoader.getDimensionsMatrix();
        // eculid makes it significant, use volume then!
        ABAMSDataMatrix dimsABAMS = new ABAMSDataMatrix( dimsMatrix, "Dimensions", new BoxDiffAdjacency() );

        // DoubleMatrix<String, String> dimsMatrix = spaceLoader.getDimensionsMatrix();
        // // eculid makes it significant, use volume then!
        // ABAMSDataMatrix dimsABAMS = new ABAMSDataMatrix( dimsMatrix, "Dimensions", new BoxDiffAdjacency() );

        StructureCatalogLoader loader = new StructureCatalogLoader();
        DoubleMatrix<String, String> nomenclatureMatrix = loader.getNomenclatureMatrix();
        ABAMSDataMatrix nomenclatureABAMS = new ABAMSDataMatrix( nomenclatureMatrix, "Nomenclature",
                new DotProductAdjacency() );

        AllenMatrixPair matrixPair;

        // get colls from ABAMS matrices?
        boolean squareMatrix = false;
        Set<String> colNames;// = new HashSet<String>( ExpressionMatrixPairFactory.getUsedCols( squareMatrix,
        // direction ) );
        colNames = new HashSet<String>( regions );

        log.info( "_____________________________" );

        // exp to space
        // use Allen matrix pair factory

        log.info( "_____________________________" );
        // exp to nomen
        matrixPair = new AllenMatrixPair( energyABAMS, nomenclatureABAMS, colNames );
        matrixPair.run();
        log.info( "_____________________________" );

        // space to nomen
        matrixPair = new AllenMatrixPair( spaceABAMS, nomenclatureABAMS, colNames );
        matrixPair.run();

        log.info( "_____________________________" );

        // // dim to nomen
        // matrixPair = new AllenMatrixPair( dimsABAMS, nomenclatureABAMS, colNames );
        // matrixPair.run();
        // log.info( "_____________________________" );
        //
        // // dim to exp
        // matrixPair = new AllenMatrixPair( dimsABAMS, energyABAMS, colNames );
        // matrixPair.run();
        // log.info( "_____________________________" );
        //
        // // dim to space
        // matrixPair = new AllenMatrixPair( dimsABAMS, spaceABAMS, colNames );
        // matrixPair.run();
        //
        // log.info( "_____________________________" );
        // // exp to vol
        // matrixPair = new AllenMatrixPair( energyABAMS, volABAMS, colNames );
        // matrixPair.run();
        // log.info( "_____________________________" );
        //
        // // space to vol
        // matrixPair = new AllenMatrixPair( spaceABAMS, volABAMS, colNames );
        // matrixPair.run();
        //
        // log.info( "_____________________________" );
        //
        // // dim to vol
        // matrixPair = new AllenMatrixPair( dimsABAMS, volABAMS, colNames );
        // matrixPair.run();
        // log.info( "_____________________________" );

    }
}
