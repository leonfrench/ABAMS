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
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.lang.time.StopWatch;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.ABAMSDataMatrix;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes;
import ubic.BAMSandAllen.SetupParameters;
import ubic.BAMSandAllen.Util;
import ubic.BAMSandAllen.AllenDataLoaders.AllenAtlasAnnotationLoader;
import ubic.BAMSandAllen.AllenDataLoaders.AllenCatalogMatrices2;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes.Direction;
import ubic.BAMSandAllen.ClassSelectors.BrainRegionClassSelector;
import ubic.BAMSandAllen.MatrixPairs.AllenMatrixPair;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenExpressionMatrixPair;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenPartialExpressionMatrixPair;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenSpacePair;
import ubic.BAMSandAllen.MatrixPairs.ExpressionMatrixPairFactory;
import ubic.BAMSandAllen.MatrixPairs.MatrixPair;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenPartialExpressionMatrixPair.RegressMatrix;
import ubic.BAMSandAllen.adjacency.BoxDiffAdjacency;
import ubic.BAMSandAllen.adjacency.CorrelationAdjacency;
import ubic.BAMSandAllen.adjacency.EuclidAdjacency;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.util.FileTools;

// use greedymutliremover
@Deprecated
public class GreedyRemover {
    private static Log log = LogFactory.getLog( GreedyRemover.class.getName() );
    MatrixPair pair;

    public GreedyRemover( MatrixPair pair ) {
        this.pair = pair;

    }

    public void removeRows( Collection<String> rows ) {
        pair.removeMatrixBDataRows( rows );
        // int count = 0;
        // for ( String row : rows ) {
        // // log.info( "Removing:" + row + " " + pair.correWithoutMatrixBDataRow( row ) );
        // if ( ++count % 1000 == 0 ) {
        // log.info( count );
        // log.info( "Removing:" + row + " " + pair.correWithoutMatrixBDataRow( row ) );
        // }
        // pair.removeMatrixBDataRowFast( row );
        // }
    }

    public void run( int iterations ) throws Exception {
        run( iterations, false );
    }

    /**
     * Iteratively remove rows from the B data matrix of the matrix pair, each time increasing the correlation the
     * maximum possible
     * 
     * @param iterations number of iterations to perform
     * @param slow indicates whether to re-compute regressions for every gene removal test
     * @throws Exception
     */
    public void run( int iterations, boolean slow ) throws Exception {
        long startTime = System.currentTimeMillis();
        StopWatch watch = new StopWatch();
        watch.start();
        StopWatch smallWatch = new StopWatch();
        int randomChoices = 0;

        for ( int i = 0; i < iterations; i++ ) {
            smallWatch.reset();
            smallWatch.start();
            // find the row that increases the correlation the most

            double baseCorrelation = pair.getCorrelation( true );
            log.info( "Base correlation:" + baseCorrelation + " size:" + pair.getMatrixBDataRows().size() );
            // pair.printDimensions();

            double bestIncrease = Double.MAX_VALUE * -1;
            String bestRow = null;
            int count = 0;
            if ( !slow && pair instanceof ConnectivityAndAllenPartialExpressionMatrixPair ) {
                ( ( ConnectivityAndAllenPartialExpressionMatrixPair ) pair ).setRegressionComputeMatrixB( false );
            }
            List<String> rows = pair.getMatrixBDataRows();
            for ( String rowName : pair.getMatrixBDataRows() ) {
                // if ( ++count % 5000 == 0 ) log.info( "Data row:" + count + " time:" + watch.getTime() );

                double withoutCorrelation = pair.correWithoutMatrixBDataRow( rowName );
                double diff = withoutCorrelation - baseCorrelation;
                if ( baseCorrelation < 0 ) {
                    diff *= -1;
                }
                if ( diff > bestIncrease ) {
                    bestRow = rowName;
                    bestIncrease = diff;
                }
            }

            if ( !slow && pair instanceof ConnectivityAndAllenPartialExpressionMatrixPair ) {
                ( ( ConnectivityAndAllenPartialExpressionMatrixPair ) pair ).setRegressionComputeMatrixB( true );
            }

            // now remove it
            // pair.removeDataRow( bestRow ); //old way
            if ( bestRow == null ) {
                bestRow = pair.getMatrixBDataRows().get( ( int ) Math.random() * rows.size() );
                log.info( "No best row found " + pair.getCorrelation( true ) );
                log.info( "Removing any gene results in correlation reduction " + pair.getCorrelation( true ) );
                break;
                // log.info( "Using random row:" + bestRow );
                // randomChoices++;
            }
            pair.removeMatrixBDataRowFast( bestRow );

            FileTools.stringToFile( bestRow + "\n", new File( SetupParameters.getDataFolder() + "LOOGenesInOrder."
                    + startTime + ".txt" ), true );
            if ( i % 1000 == 0 ) {
                // write correlation
                FileTools.stringToFile( i + "," + baseCorrelation + "\n", new File( SetupParameters.getDataFolder()
                        + "LOOResults." + startTime + ".txt" ), true );
            }
            log.info( bestRow + " changes correlation by " + bestIncrease + " time:" + ( smallWatch.getTime() / 1000 )
                    + "s total:" + ( watch.getTime() / 1000 ) + "s" );

        }
        log.info( "Random choices: " + randomChoices );
        log.info( "Start time:" + startTime );
    }

    public void checkNulls() {
        for ( String rowName : pair.getMatrixBDataRows() ) {
            if ( rowName == null ) {
                log.info( "Null name" );
            }
        }

    }

    public static void connectivityPartial( Direction direction, boolean slow, RegressMatrix regressType )
            throws Exception {
        boolean squareMatrix = false;
        if ( regressType == RegressMatrix.CONNECTIVITY ) {
            slow = false;
            log.info( "Setting to fast mode for regression on connectivity" );
        }
        // this pair is used to convert the space pair region names to BAMS names
        ConnectivityAndAllenSpacePair spacePair = new ConnectivityAndAllenSpacePair( new BrainRegionClassSelector(),
                squareMatrix, null, direction );
        spacePair.run();

        MatrixPair forR2 = new ConnectivityAndAllenPartialExpressionMatrixPair( new BrainRegionClassSelector(), true,
                false, squareMatrix, Double.NaN, "NewEnergies", direction, spacePair.getMatrixB(), regressType );

        // below is testing code
        long t;

        forR2.run();
        log.info( forR2.getCorrelation() );

        GreedyRemover remover = new GreedyRemover( forR2 );

        remover.run( forR2.getMatrixBDataRows().size() - 1, slow );

    }

    public static void dimensions( Direction direction ) throws Exception {
        AllenAtlasAnnotationLoader spaceLoader = new AllenAtlasAnnotationLoader();
        DoubleMatrix<String, String> dimsMatrix = spaceLoader.getDimensionsMatrix();
        ABAMSDataMatrix dimsABAMS = new ABAMSDataMatrix( dimsMatrix, "Dimensions", new BoxDiffAdjacency() );

        AllenCatalogMatrices2 allenMatrices = new AllenCatalogMatrices2();
        DoubleMatrix<String, String> energyMatrix = allenMatrices.getFromDisk( "NewEnergies" );
        energyMatrix = Util.logMatrix( energyMatrix, Double.NaN );
        ABAMSDataMatrix energyABAMS = new ABAMSDataMatrix( energyMatrix, "NewEnergies", new CorrelationAdjacency(
                energyMatrix ) );

        boolean squareMatrix = false;
        Set<String> colNames = new HashSet<String>( ExpressionMatrixPairFactory.getUsedCols( squareMatrix,
                direction ) );

        log.info( "_____________________________" );

        AllenMatrixPair matrixPair;

        // exp to space
        matrixPair = new AllenMatrixPair( dimsABAMS, energyABAMS, colNames );
        matrixPair.run();

        log.info( matrixPair.getCorrelation() );

        GreedyRemover remover = new GreedyRemover( matrixPair );

        // remover.run( 10 );
        remover.run( matrixPair.getMatrixBDataRows().size() - 1 );

    }

    /**
     * @param args
     */
    public static void connectivity( Direction direction, String loadFile ) throws Exception {
        // below is testing code
        long t;
        MatrixPair forR2 = new ConnectivityAndAllenExpressionMatrixPair( new BrainRegionClassSelector(), true, false,
                false, Double.NaN, "NewEnergies", direction );

        // ConnectivityAndAllenExpressionMatrixPair forR2 = new ConnectivityAndAllenExpressionMatrixPair(
        // new BrainRegionClassSelector(), true, false, false, Double.NaN, "NewEnergies" );
        forR2.run();
        log.info( forR2.getCorrelation() );

        GreedyRemover remover = new GreedyRemover( forR2 );
        if ( loadFile != null ) remover.removeRows( FileTools.getLines( loadFile ) );

        remover.run( forR2.getMatrixBDataRows().size() - 1 );

    }

    public static void space( String loadFile ) throws Exception {
        AllenAtlasAnnotationLoader spaceLoader = new AllenAtlasAnnotationLoader();
        DoubleMatrix<String, String> spaceMatrix = spaceLoader.getCenterMatrix();
        ABAMSDataMatrix spaceABAMS = new ABAMSDataMatrix( spaceMatrix, "Space", new EuclidAdjacency() );

        AllenCatalogMatrices2 allenMatrices = new AllenCatalogMatrices2();
        DoubleMatrix<String, String> energyMatrix = allenMatrices.getFromDisk( "NewEnergies" );
        energyMatrix = Util.logMatrix( energyMatrix, Double.NaN );
        ABAMSDataMatrix energyABAMS = new ABAMSDataMatrix( energyMatrix, "NewEnergies", new CorrelationAdjacency(
                energyMatrix ) );

        boolean squareMatrix = false;
        Set<String> colNames = new HashSet<String>( ExpressionMatrixPairFactory.getUsedCols( squareMatrix,
                AnalyzeBAMSandAllenGenes.Direction.OUTGOING ) );

        log.info( "_____________________________" );

        AllenMatrixPair matrixPair;

        // exp to space
        matrixPair = new AllenMatrixPair( spaceABAMS, energyABAMS, colNames );
        matrixPair.run();

        log.info( matrixPair.getCorrelation() );

        GreedyRemover remover = new GreedyRemover( matrixPair );
        remover.checkNulls();

        if ( loadFile != null ) remover.removeRows( FileTools.getLines( loadFile ) );

        // remover.run( 10 );
        remover.run( matrixPair.getMatrixBDataRows().size() - 1 );

    }

    public static void main( String[] args ) throws Exception {
        // String file = SetupParameters.getDataFolder() + "LOOGenesInOrder.Incoming.txt";
        Direction direction;
        direction = AnalyzeBAMSandAllenGenes.Direction.INCOMING;
        boolean slow = false;
        boolean regressOnConnectivity = false;
        // connectivityPartial( direction, slow, regressOnConnectivity );
        // SetupParameters.getDataFolder() + "LOOGenesInOrder.Incoming.txt" )
        // connectivity( direction, null );
        dimensions( direction );
        // space();
    }
}
