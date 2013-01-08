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

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.apache.commons.lang.time.StopWatch;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import sun.security.krb5.Config;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes;
import ubic.BAMSandAllen.RegressionVector;
import ubic.BAMSandAllen.SetupParameters;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes.Direction;
import ubic.BAMSandAllen.ClassSelectors.BrainRegionClassSelector;
import ubic.BAMSandAllen.FocusedAnalysis.TryLiteratureConnectivityMatrix;
import ubic.BAMSandAllen.MatrixPairs.AllenMatrixPair;
import ubic.BAMSandAllen.MatrixPairs.AllenMatrixPairFactory;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenExpressionMatrixPair;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenPartialExpressionMatrixPair;
import ubic.BAMSandAllen.MatrixPairs.ExpressionMatrixPairFactory;
import ubic.BAMSandAllen.MatrixPairs.MatrixPair;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenPartialExpressionMatrixPair.RegressMatrix;
import ubic.BAMSandAllen.geneFilters.GeneFilter;
import ubic.BAMSandAllen.geneFilters.NonExpFilter;
import ubic.basecode.util.FileTools;

public class GreedyMultiThreaded {
    private static Log log = LogFactory.getLog( GreedyMultiThreaded.class.getName() );
    MatrixPair pair;
    int threads;

    public GreedyMultiThreaded( MatrixPair pair, int threads ) {
        this.pair = pair;
        this.threads = threads;
    }

    public void removeRows( Collection<String> rows ) {
        pair.removeMatrixBDataRows( rows );
    }

    public void run( int iterations, boolean keepSign ) throws Exception {
        boolean slow = false;
        run( iterations, slow, keepSign );
    }

    /**
     * Iteratively remove rows from the B data matrix of the matrix pair, each time increasing the correlation the
     * maximum possible
     * 
     * @param iterations number of iterations to perform
     * @param slow indicates whether to re-compute regressions for every gene removal test
     * @throws Exception
     */
    public void run( int iterations, boolean slow, boolean keepSign ) throws Exception {
        StopWatch watch = new StopWatch();
        watch.start();
        StopWatch smallWatch = new StopWatch();
        long startTime = System.currentTimeMillis();

        double firstBaseLine = pair.getCorrelation( true );
        // make it more negative if it starts below zero
        boolean increase = firstBaseLine > 0;
        // force the increase on random runs with low correlation
        if ( Math.abs( firstBaseLine ) < 0.1 ) increase = true;

        // if we are going against the current correlation sign - eg go from positive correlation to negative
        if ( !keepSign ) {
            increase = !increase;
        }

        List<GreedyThreadRunner> runners = new LinkedList<GreedyThreadRunner>();

        // create the runners
        for ( int threadInd = 0; threadInd < threads; threadInd++ ) {
            MatrixPair pairCopy = ( MatrixPair ) deepCopy( pair );
            runners.add( new GreedyThreadRunner( pairCopy, slow, increase ) );
        }

        for ( int i = 0; i < iterations; i++ ) {
            smallWatch.reset();
            smallWatch.start();

            ExecutorService pool;
            pool = Executors.newFixedThreadPool( threads );
            // divide up the rows
            List<String> rows = pair.getMatrixBDataRows();
            List<Collection<String>> splits = split( rows, threads );

            double baseLine = pair.getCorrelation( true );

            log.info( "Base correlation:" + baseLine + " size:" + pair.getMatrixBDataRows().size() );

            // set the baseline and call the runners
            for ( int threadInd = 0; threadInd < threads; threadInd++ ) {
                GreedyThreadRunner runner = runners.get( threadInd );
                runner.setBaseline( baseLine );

                // set residual calculation triangles to match accross all threads, if we are doing partial regression
                if ( pair instanceof ConnectivityAndAllenPartialExpressionMatrixPair ) {
                    ConnectivityAndAllenPartialExpressionMatrixPair partialPair = ( ( ConnectivityAndAllenPartialExpressionMatrixPair ) pair );
                    // only do it if we are regressing on expression
                    if ( partialPair.getRegressType() == ConnectivityAndAllenPartialExpressionMatrixPair.RegressMatrix.BOTH
                            && partialPair.getRegressType() == ConnectivityAndAllenPartialExpressionMatrixPair.RegressMatrix.EXPRESSION ) {
                        RegressionVector triangles = partialPair.getTrianglesMatrixB();
                        runner.setTrianglesMatrixB( triangles );
                    }
                }

                runner.setRowsToTest( splits.get( threadInd ) );
                // log.info( "Split size:" + splits.get( threadInd ).size() );
                pool.execute( runner );
            }

            // log.info( "Waiting for threads to finish" );
            pool.shutdown();
            pool.awaitTermination( 15, TimeUnit.MINUTES );

            // go through all the runners and get results
            double bestIncrease = Double.MAX_VALUE * -1;
            String bestRow = null;
            for ( GreedyThreadRunner runner : runners ) {
                double diff = runner.getBestIncrease();
                if ( !increase ) {
                    diff *= -1;
                }
                if ( diff > bestIncrease ) {
                    bestRow = runner.getBestRow();
                    bestIncrease = diff;
                }
            }

            if ( bestRow == null ) {
                log.info( "No best row found " + pair.getCorrelation( true ) );
                log.info( "Putting remaining " + rows.size() + " genes at end of file" );
                for ( String row : rows ) {
                    outputRowToFile( startTime, row );
                    FileTools.stringToFile( i + "," + baseLine + "," + row + "\n", new File( SetupParameters
                            .getDataFolder()
                            + "LOOResults." + startTime + ".txt" ), true );
                }
                break;
            }

            pair.removeMatrixBDataRowFast( bestRow );
            for ( GreedyThreadRunner runner : runners ) {
                runner.removeRow( bestRow );
            }

            outputRowToFile( startTime, bestRow );
            // write correlation

            FileTools.stringToFile( i + "," + baseLine + "," + bestRow + "\n", new File( SetupParameters
                    .getDataFolder()
                    + "LOOResults." + startTime + ".txt" ), true );

            int eta = ( int ) ( smallWatch.getTime() / 1000 ) / 2 * pair.getMatrixBDataRows().size() / 3600;
            log.info( bestRow + " changes correlation by " + bestIncrease + " time:" + ( smallWatch.getTime() / 1000 )
                    + "s total:" + ( watch.getTime() / 1000 ) + "s estimated hours remaining:" + eta );
        }
        log.info( "Start time:" + startTime );
        FileTools.stringToFile( startTime + "\n", new File( SetupParameters.getDataFolder() + "Link." + startTime
                + ".txt" ), true );

        // make a link between it's output files and starttime - akward hack
        File jarFile = new File( this.getClass().getProtectionDomain().getCodeSource().getLocation().toURI() );
        FileTools.stringToFile( jarFile.toString() + "\n", new File( SetupParameters.getDataFolder() + "Link."
                + startTime + ".txt" ), true );

    }

    private void outputRowToFile( long startTime, String row ) throws Exception {
        FileTools.stringToFile( row + "\n", new File( SetupParameters.getDataFolder() + "LOOGenesInOrder." + startTime
                + ".txt" ), true );
    }

    /**
     * @param args
     */

    public static void connectivityLiteraturePartial( String filename, int threads, RegressMatrix regressType )
            throws Exception {
        ConnectivityAndAllenPartialExpressionMatrixPair pair = TryLiteratureConnectivityMatrix
                .getPartialLiteratureExpressionPair( filename, regressType );
        GreedyMultiThreaded remover = new GreedyMultiThreaded( pair, threads );
        boolean slow = false;
        boolean keepSign = true;
        remover.run( pair.getMatrixBDataRows().size() - 1, slow, keepSign );
    }

    public static void connectivityPartial( Direction direction, boolean slow, RegressMatrix regressType, int threads,
            boolean useVirtual, boolean removeNonExp, boolean logDistance, boolean keepSign ) throws Exception {

        MatrixPair partialMantelPair = ExpressionMatrixPairFactory.connectivityPartial( direction, slow, regressType,
                useVirtual, removeNonExp, logDistance );

        GreedyMultiThreaded remover = new GreedyMultiThreaded( partialMantelPair, threads );

        remover.run( partialMantelPair.getMatrixBDataRows().size() - 1, slow, keepSign );

    }

    public static void connectivityDirect( int threads, boolean removeNonExp ) throws Exception {

        double zeroReplace = Double.NaN;
        boolean log1p = false;
        boolean square = true;
        boolean useVirtual = false; // use virtual must be false unless the ABAMS matrix rows are made virtual too
        // only bidirectional can be used, as its symetric
        ConnectivityAndAllenExpressionMatrixPair pair = new ConnectivityAndAllenExpressionMatrixPair(
                new BrainRegionClassSelector(), true, log1p, square, zeroReplace, "NewEnergies", Direction.ANYDIRECTION );

        if ( useVirtual ) {
            pair.makeVirtualRegions();
        }

        if ( removeNonExp ) {
            pair.applyGeneFilter( new NonExpFilter() );
        }

        // pair.applyGeneFilter( new EstrogenGeneFilter() );

        pair.run();

        // default filters
        for ( GeneFilter filter : ExpressionMatrixPairFactory.getDefaultExpressionFilters() ) {
            pair.applyGeneFilter( filter );
        }
        log.info( pair.getCorrelation() );
        log.info( pair.test( 1000, false ) );

        GreedyMultiThreaded remover = new GreedyMultiThreaded( pair, threads );

        boolean keepSign = true;
        remover.run( pair.getMatrixBDataRows().size() - 1, keepSign );

    }

    public static void connectivityShuffle( Direction direction, int threads, boolean useVirtual, boolean removeNonExp,
            boolean keepSign, int seed ) throws Exception {
        ConnectivityAndAllenExpressionMatrixPair pair = ExpressionMatrixPairFactory.connectivityAndExpression(
                direction, useVirtual, removeNonExp );

        pair.shuffleDataCols( seed );
        log.info( pair.getCorrelation() );

        GreedyMultiThreaded remover = new GreedyMultiThreaded( pair, threads );

        remover.run( pair.getMatrixBDataRows().size(), keepSign );
    }

    public static void connectivityPartialShuffle( Direction direction, boolean slow, int threads, boolean useVirtual,
            boolean removeNonExp, boolean logDistance, boolean keepSign, int seed ) throws Exception {
        RegressMatrix regressType = RegressMatrix.CONNECTIVITY;
        ConnectivityAndAllenExpressionMatrixPair pair = ExpressionMatrixPairFactory.connectivityPartial( direction,
                slow, regressType, useVirtual, removeNonExp, logDistance );

        // pair.shuffleConnectivityCols( seed );
        pair.shuffleDataCols( seed );
        log.info( pair.getCorrelation() );

        GreedyMultiThreaded remover = new GreedyMultiThreaded( pair, threads );

        remover.run( pair.getMatrixBDataRows().size(), keepSign );
    }

    public static void connectivity( Direction direction, int threads, boolean useVirtual, boolean removeNonExp,
            boolean keepSign ) throws Exception {
        ConnectivityAndAllenExpressionMatrixPair pair = ExpressionMatrixPairFactory.connectivityAndExpression(
                direction, useVirtual, removeNonExp );

        GreedyMultiThreaded remover = new GreedyMultiThreaded( pair, threads );

        remover.run( pair.getMatrixBDataRows().size(), keepSign );
    }

    List<Collection<String>> split( List<String> list, int ways ) {
        int size = list.size() / ways;
        List<Collection<String>> result = new LinkedList<Collection<String>>();
        for ( int i = 0; i < ways; i++ ) {
            if ( i == ways - 1 ) { // its the last one
                result.add( list.subList( i * size, list.size() ) );
            } else {
                result.add( list.subList( i * size, ( i + 1 ) * size ) );
            }
        }
        return result;
    }

    public static void volumes( Direction direction, int threads, boolean keepSign, boolean removeNonExp )
            throws Exception {
        AllenMatrixPair matrixPair = AllenMatrixPairFactory.getVolumeAndExpressionPair( direction, removeNonExp );

        GreedyMultiThreaded remover = new GreedyMultiThreaded( matrixPair, threads );

        remover.run( matrixPair.getMatrixBDataRows().size() - 1, keepSign );
    }

    public static void space( Direction direction, int threads, boolean removeNonExp, boolean logDistance,
            boolean keepSign ) throws Exception {

        AllenMatrixPair matrixPair = AllenMatrixPairFactory.getSpaceAndExpressionPair( direction, removeNonExp,
                logDistance );

        GreedyMultiThreaded remover = new GreedyMultiThreaded( matrixPair, threads );

        remover.run( matrixPair.getMatrixBDataRows().size(), keepSign );
    }

    public static void volume( Direction direction, int threads, boolean removeNonExp ) throws Exception {
        AllenMatrixPair matrixPair = AllenMatrixPairFactory.getVolumeAndExpressionPair( direction, removeNonExp );

        GreedyMultiThreaded remover = new GreedyMultiThreaded( matrixPair, threads );

        boolean keepSign = true;
        remover.run( matrixPair.getMatrixBDataRows().size(), keepSign );
    }

    public static void removeRegions( MatrixPair pair, int threads ) throws Exception {
        pair.switchMatrices();
        GreedyMultiThreaded remover = new GreedyMultiThreaded( pair, threads );
        boolean keepSign = true;
        remover.run( pair.getMatrixBDataRows().size(), keepSign );

    }

    public static void main( String[] args ) throws Exception {
        Direction direction;
        direction = AnalyzeBAMSandAllenGenes.Direction.INCOMING;
        boolean slow = false;
        boolean logDistance = true;
        RegressMatrix regressType = RegressMatrix.CONNECTIVITY;
        int threads = 8;
        boolean virtual = true;
        boolean removeNonExp = true;
        boolean keepSign = true;

        // nice java -server -classpath
        // icu4j-3.4.jar:xercesImpl-2.8.1.jar:BAMSandAllen_fat.literatureConnect.expression.jar -Xmx27000M
        // ubic.BAMSandAllen.optimize.GreedyMultiThreaded

        // ExpressionMatrixPairFactory.addFilter( new ABA1000ConflictsGeneFilter() );

        String base = SetupParameters.config.getString( "whitetext.iteractions.matricesFolder" );
        String filename = base;
        filename += "Positives.WhiteTextUnseen.matrix.txt.propigated"; // 115 columns used

        connectivityLiteraturePartial( filename, threads, regressType );

        // connectivityPartial( direction, slow, regressType, threads, virtual, removeNonExp, logDistance, keepSign );

        //
        // space( direction, threads, removeNonExp, logDistance, keepSign );
        // volume( direction, threads, removeNonExp );
        // connectivity( direction, threads, virtual, removeNonExp, keepSign );
        // connectivityDirect( threads, removeNonExp );

        // run with incoming axonongenisis GO groups
        // direction = AnalyzeBAMSandAllenGenes.Direction.INCOMING;
        // connectivityPartialShuffle( direction, slow, threads, virtual, removeNonExp, logDistance, keepSign, 2 );

        // for ( int seed = 0; seed < 10; seed++ ) {
        // connectivityPartialShuffle( direction, slow, threads, virtual, removeNonExp, logDistance, keepSign, seed );
        // // connectivityShuffle( direction, threads, virtual, removeNonExp, keepSign, seed );
        // }

    }

    public static Object deepCopy( Object oldObj ) throws Exception {
        ByteArrayOutputStream bos = new ByteArrayOutputStream();
        ObjectOutputStream oos = new ObjectOutputStream( bos );
        oos.writeObject( oldObj );
        oos.flush();
        ObjectInputStream ois = new ObjectInputStream( new ByteArrayInputStream( bos.toByteArray() ) );
        return ois.readObject();
    }
}
