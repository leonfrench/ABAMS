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

import java.io.File;
import java.io.Serializable;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.ABAMSDataMatrix;
import ubic.BAMSandAllen.SetupParameters;
import ubic.BAMSandAllen.Util;
import ubic.BAMSandAllen.geneFilters.GeneFilter;
import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.math.CorrelationStats;
import ubic.basecode.math.DescriptiveWithMissing;
import ubic.basecode.util.FileTools;
import cern.colt.list.DoubleArrayList;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.jet.math.Functions;

/*
 * two pairs of matrices, same rows but different columns
 * 
 * 
 * the columns represent brain regions
 */
public abstract class MatrixPair implements Serializable {
    private static Log log = LogFactory.getLog( MatrixPair.class.getName() );
    String baseName;

    boolean isInSameSpace;

    boolean spearman = false;

    // DoubleMatrix<String, String> matrixA;
    // DoubleMatrix<String, String> matrixB;
    ABAMSDataMatrix matrixA;
    ABAMSDataMatrix matrixB;

    public MatrixPair() {
        baseName = SetupParameters.getDataFolder() + this.getClass().getSimpleName();
        isInSameSpace = false;
    }

    public String getName() {
        return this.getClass().getSimpleName();
    }

    public String getBaseName() {
        return baseName;
    }

    // needs to be renamed
    public abstract double run() throws Exception;

    public abstract Set<String> convertBNametoA( String bName );

    public abstract Set<String> convertANametoB( String aName );

    public ABAMSDataMatrix getMatrixA() {
        return matrixA;
    }

    public ABAMSDataMatrix getMatrixB() {
        return matrixB;
    }

    public double[] getCorrelationTriangle( ABAMSDataMatrix matrix ) {
        DoubleMatrix<String, String> adjacency = matrix.getAdjacency();
        double[] vecTri = Util.getTriangle( adjacency );
        return vecTri;
    }

    public void applyMatrixBRowFilter( GeneFilter filter ) throws Exception {
        ABAMSDataMatrix matrix = getMatrixB();
        matrixB = matrix.applyRowFilter( filter );
        // String originalName= matrixB.getName();
        // originalName+="."+ filter.getName();
        // matrixB.setName( originalName );
    }

    /*
     * ensures the col names are the same order and same entries, matrixB col names are converted to matrixA names
     */
    public void sameSpace() {
        DoubleMatrix<String, String> newMatrixB = new DenseDoubleMatrix<String, String>( matrixB.rows(), matrixA
                .columns() );
        newMatrixB.setColumnNames( matrixA.getColNames() );
        newMatrixB.setRowNames( matrixB.getRowNames() );

        for ( String matrixAColumn : matrixA.getColNames() ) {
            Set<String> matrixBColumns = convertANametoB( matrixAColumn );
            // we may still have been mapped to a brain region with no data
            // requires mapping and matrix data!
            matrixBColumns.retainAll( matrixB.getColNames() );

            if ( matrixBColumns.size() > 1 ) log.info( "Merging " + matrixBColumns + " into " + matrixAColumn );

            // merge them all in - all of these B columns gota be put into a single spot in the new Matrix
            for ( String matrixBColumn : matrixBColumns ) {
                // do the averaging into the BAMSRegion
                int oldIndex = matrixB.getColIndexByName( matrixBColumn );
                int newIndex = newMatrixB.getColIndexByName( matrixAColumn );

                // the column of data to be merged
                double[] matrixColumn = matrixB.getColumn( oldIndex );

                for ( int row = 0; row < matrixColumn.length; row++ ) {
                    double current = newMatrixB.get( row, newIndex );
                    // average things out
                    double value = current + matrixColumn[row] / ( ( double ) matrixBColumns.size() );
                    newMatrixB.set( row, newIndex, value );
                }
            }
        }
        matrixB = matrixB.replaceMatrix( newMatrixB );
        isInSameSpace = true;
    }

    public void testZeroes() {
        log.info( matrixA.getName() + " Zeroes columns:" + Util.findZeroColumns( matrixA ).size() );
        log.info( matrixB.getName() + " Zeroes columns:" + Util.findZeroColumns( matrixB ).size() );
        log.info( matrixA.getName() + " Zero rows:" + Util.findZeroRows( matrixA ).size() );
        log.info( matrixB.getName() + " Zeroes rows:" + Util.findZeroRows( matrixB ).size() );
    }

    public void testZeroesVerbose() {
        log.info( matrixA.getName() + " Zeroes:" + Util.findZeroColumns( matrixA ) );
        log.info( matrixB.getName() + " Zeroes:" + Util.findZeroColumns( matrixB ) );
    }

    /*
     * beware
     */
    public void switchMatrices() {
        log.warn( "SWitching matrices" );
        ABAMSDataMatrix temp = matrixA;
        matrixA = matrixB;
        matrixB = temp;
    }

    public double testSwitch( int iterations, boolean quiet ) throws Exception {
        switchMatrices();
        double result = test( iterations, quiet );
        // switch it back
        switchMatrices();
        return result;
    }

    public void testBoth( int iterations, boolean verbose ) throws Exception {
        test( iterations, verbose );
        testSwitch( iterations, verbose );
    }

    public double test( int iterations ) throws Exception {
        return test( iterations, true );
    }

    public double test( int iterations, boolean verbose ) throws Exception {
        DoubleMatrix<String, String> aCorrelations = matrixA.getAdjacency(); // Util.correlateColumns( matrixA,
        // spearman );
        DoubleMatrix<String, String> bCorrelations = matrixB.getAdjacency(); // Util.correlateColumns( matrixB,
        // spearman );
        // double real = Math.abs( CorrelationStats.correl( getCorrelationTriangle( matrixB ),
        // getCorrelationTriangle( matrixA ) ) );
        double real = Math.abs( getCorrelation( true ) );
        if ( verbose ) log.info( "Using " + real + " as threshold" );
        double[] aCorrelationTriangle = matrixA.getTriangle();
        int lowerScore = 0;
        // shuffle be
        List<String> names = new LinkedList<String>( bCorrelations.getColNames() );
        Random random = new Random( 1 );

        for ( int i = 0; i < iterations; i++ ) {
            if ( verbose && i % 1000 == 0 ) log.info( "Shuffle " + i );
            Collections.shuffle( names, random );
            DoubleMatrix<String, String> shuffledCor = new DenseDoubleMatrix<String, String>( names.size(), names
                    .size() );
            shuffledCor.setColumnNames( names );
            shuffledCor.setRowNames( names );

            for ( int row = 0; row < shuffledCor.rows(); row++ ) {
                for ( int column = 0; column < shuffledCor.columns(); column++ ) {
                    String colName = names.get( column );
                    String rowName = names.get( row );

                    shuffledCor.set( row, column, bCorrelations.getByKeys( rowName, colName ) );
                }
            }

            double shuffledResult = CorrelationStats.correl( Util.getTriangle( shuffledCor ), aCorrelationTriangle );

            if ( Math.abs( shuffledResult ) >= real ) {
                if ( verbose ) log.info( shuffledResult );
                lowerScore++;
            }
        }

        log.info( lowerScore + " of " + iterations + " random sets have correlation with absolute value above " + real );
        double result = lowerScore / ( double ) iterations;
        if ( lowerScore == 0 ) result = 0;
        return result;
    }

    public void printMappingRelations() {
        int aToManyB = 0;
        int bToManyA = 0;
        int noMappingA = 0;
        int noMappingB = 0;

        for ( String aName : matrixA.getColNames() ) {
            Set<String> bNames = convertANametoB( aName );
            if ( bNames != null ) {
                if ( bNames.size() > 1 ) {
                    aToManyB++;
                    log.info( matrixA.getName() + " " + aName + " -> " + bNames );
                }
            }
            if ( bNames == null || bNames.isEmpty() ) {
                // log.info( matrixA.getName() + " " + aName + " has no mapping" );
                noMappingA++;
            }
        }

        // bName -> aName
        for ( String bName : matrixB.getColNames() ) {
            Set<String> aNames = convertBNametoA( bName );
            if ( aNames != null ) {
                if ( aNames.size() > 1 ) {
                    bToManyA++;
                    log.info( matrixB.getName() + " " + bName + " -> " + aNames );
                }
            }
            if ( aNames == null || aNames.isEmpty() ) {
                // log.info( matrixB.getName() + " " + bName + " has no mapping" );
                noMappingB++;
            }
        }
        log.info( matrixB.getName() + " columns mapped to many " + matrixA.getName() + ": " + bToManyA );
        log.info( matrixA.getName() + " columns mapped to many " + matrixB.getName() + ": " + aToManyB );
        log.info( matrixB.getName() + " columns had no mapping " + matrixA.getName() + ": " + noMappingB );
        log.info( matrixA.getName() + " columns had no mapping " + matrixB.getName() + ": " + noMappingA );
    }

    public void slimMatrices() {
        // run it twice to get both matrices slimmed
        slimMatricesOnce();
        slimMatricesOnce();
    }

    /*
     * reduce the matrix to the columns that survive the mapping
     */
    public final void slimMatricesOnce() {
        log.info( "Before" );
        printDimensions();
        // get seen, then retain those with current matrix!!!
        // aName -> bName
        Set<String> bMapped = new HashSet<String>();
        for ( String aName : matrixA.getColNames() ) {
            Set<String> bNames = convertANametoB( aName );
            if ( bNames != null ) {
                bMapped.addAll( convertANametoB( aName ) );
            }
        }

        // bName -> aName
        Set<String> aMapped = new HashSet<String>();
        for ( String bName : matrixB.getColNames() ) {
            Set<String> aNames = convertBNametoA( bName );
            if ( aNames != null ) {
                aMapped.addAll( aNames );
            }
        }

        // we may have been mapped to a name that has no data
        bMapped.retainAll( matrixB.getColNames() );
        aMapped.retainAll( matrixA.getColNames() );

        // we only realy need to slim A
        matrixA = matrixA.retainColumns( aMapped );
        matrixB = matrixB.retainColumns( bMapped );

        log.info( "After" );

        log.info( "Removing all zero rows" );

        // if this is done, it will increase the weight of up-propigated connections
        // matrixA = matrixA.removeRows( Util.findZeroRows( matrixA ) );

        matrixB = matrixB.removeRows( Util.findZeroRows( matrixB ) );

        printDimensions();
    }

    public boolean hasNaNCorrelation( Collection<String> removeRows ) {
        Set<String> removeRowsSet = new HashSet<String>( removeRows );
        ABAMSDataMatrix matrixBforCor = matrixB.removeRows( removeRowsSet );
        // log.info( matrixBName + " = " + matrixBforCor.rows() + " x " + matrixBforCor.columns() );
        DoubleMatrix<String, String> bCorrelations = Util.correlateColumns( matrixBforCor, spearman );
        return Util.countValues( bCorrelations, Double.NaN ) > 0;
    }

    public DoubleMatrix<String, String> getAdjacencyB() {
        return matrixB.getAdjacency();
    }

    public DoubleMatrix<String, String> getAdjacencyA() {
        return matrixA.getAdjacency();
    }

    public double getCorrelation() throws Exception {
        return getCorrelation( matrixB, false );
    }

    public double getCorrelation( boolean fast ) throws Exception {
        return getCorrelation( matrixB, fast );
    }

    public double getCorrelation( ABAMSDataMatrix matrixBforCor, boolean fast ) throws Exception {
        String matrixAName = matrixA.getName();
        String matrixBName = matrixB.getName();
        if ( !matrixA.getColNames().equals( matrixB.getColNames() ) ) {
            throw new Exception( "Col names do not match" );
        }

        if ( !fast ) {
            // two square matrices - size is based on BAMS mapped terms but Allen could be used (how to merge connection
            // profiles?)
            // don't use spearman on connections - its binary
            DoubleMatrix<String, String> aCorrelations = matrixA.getAdjacency();
            // very slow at this call
            DoubleMatrix<String, String> bCorrelations = matrixBforCor.getAdjacency();

            // min and max code should be refactored
            DoubleFactory2D fact = DoubleFactory2D.dense;
            DoubleMatrix2D a = fact.make( aCorrelations.getRawMatrix() );
            DoubleMatrix2D b = fact.make( bCorrelations.getRawMatrix() );
            // set diagonals to first non diagonal cell
            for ( int i = 0; i < a.rows(); i++ ) {
                a.set( i, i, a.get( 0, 1 ) );
                b.set( i, i, b.get( 0, 1 ) );
            }
            // diagonal is zero so it gives the wrong min, fix later
            log.info( matrixAName + " Max:" + a.aggregate( Functions.max, Functions.identity ) );
            log.info( matrixAName + " Min:" + a.aggregate( Functions.min, Functions.identity ) );
            log.info( matrixBName + " Max:" + b.aggregate( Functions.max, Functions.identity ) );
            log.info( matrixBName + " Min:" + b.aggregate( Functions.min, Functions.identity ) );

            try {
                Util.writeImage( baseName + "." + matrixAName + ".Correlation.png", aCorrelations );
                Util.writeImage( baseName + "." + matrixBName + ".Correlation.png", bCorrelations );

                DoubleMatrix<String, String> forPlot = Util.getTrianglePairedPoints( aCorrelations, bCorrelations,
                        matrixAName, matrixBName );
                Util.writeRTable( baseName + "." + matrixAName + "." + matrixBName + ".forPlotting.txt", forPlot );

                Util.writeRTable( baseName + "." + matrixAName + ".Correlation.txt", aCorrelations );
                Util.writeRTable( baseName + "." + matrixBName + ".Correlation.txt", bCorrelations );

            } catch ( Exception e ) {
                e.printStackTrace();
            }
        }

        double[] aVecTri = matrixA.getTriangle();
        double[] bVecTri = matrixBforCor.getTriangle();

        if ( !fast ) {
            log.info( matrixBName + " triangle Size:" + bVecTri.length );
            log.info( matrixAName + " triangle Size:" + aVecTri.length );
        }

        if ( bVecTri.length != aVecTri.length ) throw new RuntimeException( "Error triangles different lengths" );

        return CorrelationStats.correl( aVecTri, bVecTri );
    }

    public void printDimensions() {
        log.info( matrixA.getDimensionString() );
        log.info( matrixB.getDimensionString() );
    }

    public void writeRMatrices() throws Exception {
        Util.writeRTable( baseName + "." + matrixA.getName() + ".txt", matrixA );
        Util.writeRTable( baseName + "." + matrixB.getName() + ".txt", matrixB );
    }

    public void writeImages() throws Exception {
        Util.writeImage( baseName + "." + matrixA.getName() + ".png", matrixA );
        Util.writeImage( baseName + "." + matrixB.getName() + ".png", matrixB );
    }

    public double correWithoutMatrixBDataRow( String rowName ) {
        double[] matrixATri = matrixA.getTriangle();
        double[] matrixBTri = matrixB.getTriangleWithoutRow( rowName );
        double real = CorrelationStats.correl( matrixATri, matrixBTri );
        return real;
    }

    /*
     * operation on the matrix, should be in another class
     */
    public ABAMSDataMatrix removeZeroLines( ABAMSDataMatrix matrix ) {
        Set<String> removeCols = Util.findZeroColumns( matrix );
        log.info( matrixA.getName() + " removing zero cols and rows =" + removeCols.size() );
        matrix = matrix.removeZeroColumns();
        matrix = matrix.removeRows( removeCols );
        return matrix;
    }

    /*
     * make it so matrix A and B have the same column names, do it by averaging or whatever..
     */
    public List<String> getMatrixBDataRows() {
        return matrixB.getRowNames();
    }

    public void removeMatrixBDataRows( Collection<String> rows ) {
        matrixB = matrixB.removeRows( rows );
    }

    public void removeMatrixADataRows( Collection<String> rows ) {
        matrixA = matrixA.removeRows( rows );
    }

    public void setMatrixBDataRows( Collection<String> rows ) {
        List<String> currentRows = new LinkedList<String>( matrixB.getRowNames() );
        currentRows.removeAll( rows );
        log.info( "Removing " + currentRows.size() + " data rows" );
        removeMatrixBDataRows( currentRows );
    }

    public Map<String, String> testDataRows( List<String> rows, int testSamples, boolean degree, boolean spearman )
            throws Exception {
        // result - probably should be an object
        Map<String, String> results = new HashMap<String, String>();
        // backup matrix
        ABAMSDataMatrix originalMatrixB = matrixB;

        if ( Util.intersect( rows, matrixB.getRowNames() ).isEmpty() ) {
            throw new RuntimeException( "Error, no rows remain" );
        }
        setMatrixBDataRows( rows );

        // size
        int size = rows.size();
        results.put( "size", size + "" );
        // correlation
        String statisticString = "correlation";
        if ( degree ) statisticString = "degree correlation";
        if ( spearman && degree ) statisticString = "degree Spearman correlation";

        double statistic;
        if ( degree ) {
            statistic = getFlattenedCorrelation( spearman );
        } else {
            // if ( spearman ) throw new RuntimeException( "Spearman not supported" );
            statistic = getCorrelation();
        }
        results.put( statisticString, statistic + "" );
        // p-value shuffle
        String spearmanString = "";
        if ( spearman ) spearmanString = " Spearman";
        if ( !degree ) results.put( "p-value shuffle" + spearmanString, test( testSamples ) + "" );
        // p-value for resampling
        Random r = new Random( 1 );
        int greaterHits = 0;
        int lesserHits = 0;
        // slow
        File f = new File( SetupParameters.getDataFolder() + statisticString + ".distro." + size + ".txt" );
        f.delete();

        DoubleArrayList sampleHistory = new DoubleArrayList();

        if ( !rows.equals( getMatrixBDataRows() ) ) {
            for ( int i = 0; i < testSamples; i++ ) {
                matrixB = originalMatrixB;
                // randomly choose a set of the same size
                List<String> rownames = new LinkedList<String>( matrixB.getRowNames() );
                Collections.shuffle( rownames, r );
                setMatrixBDataRows( rownames.subList( 0, size ) );
                double sampleStatistic;

                if ( !degree )
                    sampleStatistic = getCorrelation();
                else
                    sampleStatistic = getFlattenedCorrelation( spearman );

                FileTools.stringToFile( sampleStatistic + "\n", f, true );
                sampleHistory.add( sampleStatistic );

                if ( sampleStatistic < statistic ) {
                    lesserHits++;
                    log.info( "Sample " + statisticString + ":" + sampleStatistic + " less than " + statistic );
                }
                if ( sampleStatistic > statistic ) {
                    greaterHits++;
                    log.info( "Sample " + statisticString + ":" + sampleStatistic + " greater than " + statistic );
                }
            }
            double mean = DescriptiveWithMissing.mean( sampleHistory );
            log.info( "Sample " + statisticString + " greater than " + statistic + ", " + greaterHits + " times" );
            log.info( "Sample " + statisticString + " less than " + statistic + ", " + lesserHits + " times" );
            results.put( "average resample" + spearmanString, mean + "" );
            results.put( "standard deviation resample" + spearmanString, Math.sqrt( DescriptiveWithMissing
                    .sampleVariance( sampleHistory, mean ) )
                    + "" );
            double finalHits = 0;
            boolean leftTail = mean > statistic;
            results.put( "left tail" + spearmanString, leftTail + "" );
            // Choose tail based on center
            if ( !leftTail )
                finalHits = greaterHits;
            else
                finalHits = lesserHits;
            results.put( "p-value resample" + spearmanString, finalHits / ( double ) testSamples + "" );
            // put it back to the original matrix
        }
        matrixB = originalMatrixB;
        return results;
    }

    /**
     * Sums all of matrixB values and sums all of matrixA values and does a correlation
     * 
     * @return
     */

    public double getFlattenedCorrelation( boolean spearman ) throws Exception {
        boolean removeNan = true;
        DoubleMatrix<String, String> matrixASums = Util.columnSums( matrixA, removeNan );
        DoubleMatrix<String, String> matrixBSums = Util.columnSums( matrixB, removeNan );

        // are they in the right order?
        if ( !matrixASums.getColNames().equals( matrixB.getColNames() ) )
            throw new RuntimeException( "Error column names do not match" );
        
        double[] matrixAVec = matrixASums.getRow( 0 );
        double[] matrixBVec = matrixBSums.getRow( 0 );
        if ( !spearman ) {
            return CorrelationStats.correl( matrixAVec, matrixBVec );
        } else {
            return Util.spearmanCorrel( matrixAVec, matrixBVec );
        }
    }

    public void removeMatrixBDataRowFast( String row ) {
        // List<String> rows = new LinkedList<String>();
        // rows.add( row );
        // boolean clone = false;
        matrixB = matrixB.copyWithoutRow( row );
        // matrixB = matrixB.removeRows( rows, clone );
    }

}
