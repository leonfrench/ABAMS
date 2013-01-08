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
package ubic.BAMSandAllen;

import java.io.File;
import java.io.FileWriter;
import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.AllenDataLoaders.ImageSeriesInfoLoader;
import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.graphics.MatrixDisplay;
import ubic.basecode.io.writer.MatrixWriter;
import ubic.basecode.math.CorrelationStats;
import ubic.basecode.math.Distance;
import ubic.basecode.util.FileTools;
import cern.colt.list.DoubleArrayList;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;

public class Util {
    private static Log log = LogFactory.getLog( Util.class.getName() );

    // find the all zero columns of a matrix
    public static Set<String> findZeroColumns( DoubleMatrix<String, String> input ) {
        Set<String> result = new HashSet<String>();
        outer: for ( int i = 0; i < input.columns(); i++ ) {
            double[] iCol = input.getColumn( i );
            for ( double d : iCol ) {
                if ( d != 0d ) continue outer;
            }
            result.add( input.getColName( i ) );
        }
        return result;
    }

    public static Set<String> findZeroRows( DoubleMatrix<String, String> input ) {
        Set<String> result = new HashSet<String>();
        outer: for ( int i = 0; i < input.rows(); i++ ) {
            double[] iCol = input.getRow( i );
            for ( double d : iCol ) {
                if ( d != 0d ) continue outer;
            }
            result.add( input.getRowName( i ) );
        }
        return result;
    }

    public static int countValues( DoubleMatrix<String, String> input, double match ) {
        int result = 0;
        boolean isNaN = Double.isNaN( match );
        for ( int row = 0; row < input.rows(); row++ ) {
            for ( int col = 0; col < input.columns(); col++ ) {
                if ( isNaN && Double.isNaN( input.get( row, col ) ) )
                    result++;
                else if ( input.get( row, col ) == match ) result++;
            }
        }
        return result;
    }

    public static int countNaNs( double[] array ) {
        int count = 0;
        for ( double d : array )
            if ( Double.isNaN( d ) ) count++;
        return count;
    }

    /*
     * remove rows from a matrix
     */
    public static DoubleMatrix<String, String> removeRows( DoubleMatrix<String, String> input,
            Collection<String> removeRows ) {
        // fix - needs index?
        Set<String> rowsToKeep = new HashSet<String>( input.getRowNames() );
        rowsToKeep.removeAll( removeRows );

        DoubleMatrix<String, String> output = new DenseDoubleMatrix<String, String>( rowsToKeep.size(), input.columns() );
        output.setColumnNames( input.getColNames() );
        int index = 0;
        for ( String row : input.getRowNames() ) {
            if ( removeRows.contains( row ) ) continue;
            output.setRowName( row, index++ );
            
            int newRowIndex = output.getRowIndexByName( row );
            int oldRowIndex = input.getRowIndexByName( row );
            for ( int col = 0; col < output.columns(); col++ ) {
                output.set( newRowIndex, col, input.get( oldRowIndex, col ) );
            }
        }

        return output;
    }

    // http://warrenseen.com/blog/2006/03/13/how-to-calculate-standard-deviation/
    public static double variance( double[] population ) {
        long n = 0;
        double mean = 0;
        double s = 0.0;

        for ( double x : population ) {
            n++;
            double delta = x - mean;
            mean += delta / n;
            s += delta * ( x - mean );
        }
        // if you want to calculate std deviation
        // of a sample change this to (s/(n-1))
        return ( s / n );
    }

    public static double spearmanCorrel( double[] a, double[] b ) {
        return Distance.spearmanRankCorrelation( new DoubleArrayList( a ), new DoubleArrayList( b ) );
    }

    public static DoubleMatrix<String, String> correlateColumns( DoubleMatrix<String, String> input, boolean spearman ) {
        int size = input.columns();
        DoubleMatrix<String, String> result = new DenseDoubleMatrix<String, String>( size, size );
        result.setColumnNames( input.getColNames() );
        result.setRowNames( input.getColNames() );
        for ( int i = 0; i < size; i++ ) {
            double[] iCol = input.getColumn( i );
            for ( int j = i; j < size; j++ ) {
                double[] jCol = input.getColumn( j );
                double value;
                if ( spearman ) {
                    value = spearmanCorrel( iCol, jCol );
                } else { // compute Pearson
                    value = CorrelationStats.correl( iCol, jCol );
                }

                result.set( i, j, value );
                result.set( j, i, value );
            }
        }
        return result;
    }

    // upper right
    public static DoubleMatrix<String, String> getTrianglePairedPoints( DoubleMatrix<String, String> a,
            DoubleMatrix<String, String> b, String nameA, String nameB ) throws Exception {
        if ( a.columns() != a.rows() ) throw new RuntimeException( "Not a square matrix" );
        if ( b.columns() != b.rows() ) throw new RuntimeException( "Not a square matrix" );
        if ( a.columns() != b.rows() ) throw new RuntimeException( "Not matching matrices" );

        int size = a.columns();
        DoubleMatrix<String, String> result = new DenseDoubleMatrix<String, String>( size * ( size - 1 ) / 2, 2 );
        result.addColumnName( nameA );
        result.addColumnName( nameB );

        for ( int i = 0; i < size; i++ ) {
            for ( int j = i + 1; j < size; j++ ) {
                String colName = a.getColName( i );
                String rowName = a.getColName( j );
                String newRowName = colName + "-" + rowName;
                result.addRowName( newRowName );
                result.setByKeys( newRowName, nameA, a.get( i, j ) );
                result.setByKeys( newRowName, nameB, b.get( i, j ) );
            }
        }
        return result;
    }

    // upper right
    public static double[] getTriangle( DoubleMatrix<String, String> input ) {
        if ( input.columns() != input.rows() ) throw new RuntimeException( "Not a square matrix" );
        int size = input.columns();
        double[] result = new double[size * ( size - 1 ) / 2];
        int spot = 0;
        for ( int i = 0; i < size; i++ ) {
            for ( int j = i + 1; j < size; j++ ) {
                result[spot++] = input.get( i, j );
            }
            // for ( int j = 0; j < i; j++ ) {
            // result[spot++] = input.get( i, j );
            // }
        }
        return result;
    }

    /*
     * zeroReplacement happens before log
     */
    public static DoubleMatrix<String, String> logMatrix( DoubleMatrix<String, String> input, double zeroReplacement ) {
        DoubleMatrix<String, String> result = input.copy();
        for ( int i = 0; i < result.rows(); i++ ) {
            for ( int j = 0; j < result.columns(); j++ ) {
                double value = result.get( i, j );
                if ( value == 0 ) {
                    value = zeroReplacement;
                }
                result.set( i, j, Math.log( value ) );
            }
        }
        return result;
    }

    public static DoubleMatrix<String, String> log1pMatrix( DoubleMatrix<String, String> input ) {
        DoubleMatrix<String, String> result = input.copy();
        for ( int i = 0; i < result.rows(); i++ ) {
            for ( int j = 0; j < result.columns(); j++ ) {
                double value = result.get( i, j );
                result.set( i, j, Math.log1p( value ) );
            }
        }
        return result;
    }

    public static DoubleMatrix<String, String> columnSums( DoubleMatrix<String, String> input ) {
        boolean removeNaN = false;
        return columnSums( input, removeNaN );
    }

    public static DoubleMatrix<String, String> columnSums( DoubleMatrix<String, String> input, boolean removeNaN ) {
        DoubleMatrix<String, String> result = new DenseDoubleMatrix<String, String>( 1, input.columns() );
        result.setColumnNames( input.getColNames() );
        result.addRowName( "Sums" );
        double sum;
        for ( int j = 0; j < input.columns(); j++ ) {
            double[] d = input.getColumn( j );
            sum = sum( d, removeNaN );
            result.set( 0, j, sum );
        }
        return result;
    }

    public static double sum( double[] array ) {
        return sum( array, false );
    }

    public static double sum( double[] array, boolean removeNaN ) {
        double result = 0;
        for ( double d : array ) {
            if ( removeNaN && Double.isNaN( d ) ) {
                // skip it
                continue;
            } else {
                result += d;
            }
        }

        return result;
    }

    public static void main( String args[] ) {
        DoubleMatrix<String, String> matrix = new DenseDoubleMatrix<String, String>( 3, 3 );
        matrix.set( 0, 0, -1d );
        matrix.set( 0, 1, 1d );
        matrix.set( 0, 2, 2d );
        matrix.set( 1, 0, 4d );
        matrix.set( 1, 1, -1d );
        matrix.set( 1, 2, 5d );
        matrix.set( 2, 0, 6d );
        matrix.set( 2, 1, 7d );
        matrix.set( 2, 2, -1d );
        matrix.addRowName( "a" );
        matrix.addRowName( "b" );
        matrix.addRowName( "c" );
        matrix.setColumnNames( matrix.getRowNames() );

        System.out.println( getTriangle( matrix )[0] );
        System.out.println( getTriangle( matrix )[1] );
        System.out.println( getTriangle( matrix )[2] );

        System.out.println( matrix.toString() );
        matrix = logMatrix( matrix, 0.000001 );
        System.out.println( matrix.toString() );

        System.out.println( "log -1 " + Math.log( -1 ) );
        System.out.println( "log 0 " + Math.log( 0 ) );
        System.out.println( "log NaN " + Math.log( Double.NaN ) );

    }

    public static void writeRTable( String filename, DoubleMatrix<String, String> matrix ) throws Exception {
        // write it out for R
        DoubleMatrix<String, String> writeMatrix = matrix.copy();

        // quote them for R
        List<String> newColNames = quoteList( matrix.getColNames() );
        writeMatrix.setColumnNames( newColNames );

        List<String> newRowNames = quoteList( matrix.getRowNames() );
        writeMatrix.setRowNames( newRowNames );

        FileWriter fOut = new FileWriter( filename );
        MatrixWriter matWriter = new MatrixWriter( fOut );
        matWriter.setTopLeft( "" );
        matWriter.writeMatrix( writeMatrix, true );
        fOut.close();
    }

    private static List<String> quoteList( List<String> oldNames ) {
        List<String> newColNames = new LinkedList<String>();
        for ( String colName : oldNames ) {
            newColNames.add( "\"" + colName + "\"" );
        }
        return newColNames;
    }

    public static void writeImage( String filename, DoubleMatrix<String, String> matrix ) throws Exception {
        MatrixDisplay matDisplay = new MatrixDisplay( matrix );
        matDisplay.setLabelsVisible( true );
        matDisplay.saveImage( filename );
    }

    public static Set<String> getNCBIDSFromGeneNames( String filename, List<String> geneNames ) throws Exception {
        ImageSeriesInfoLoader imageInfo = new ImageSeriesInfoLoader();
        List<String> rowNames = getRowNamesFromGeneNames( filename, geneNames );
        Set<String> result = new HashSet<String>();
        for ( String rowName : rowNames ) {
            String NCBID = imageInfo.getNCBIIDFromRowName( rowName ) + "";
            result.add( NCBID );
        }
        return result;
    }

    /**
     * Given a file with gene names, it will filter a list of rowNames for those that have the gene names as prefixes.
     * It is usefull for the Allen gene name lists of non expressors or ubiquituosly expressed genes.
     * 
     * @param filename - input file
     * @param rowNames - set of rownames to match against
     */
    public static List<String> getRowNamesFromGeneNames( String filename, List<String> rowNames ) throws Exception {
        List<String> genes = FileTools.getLines( filename );

        int noHits = 0;
        List<String> resultRows = new LinkedList<String>();
        for ( String gene : genes ) {
            boolean hit = false;
            for ( String row : rowNames ) {
                if ( row.startsWith( gene + "[" ) ) {
                    resultRows.add( row );
                    hit = true;
                }
            }
            if ( !hit ) noHits++;
        }
        log.info( "Genes from file:" + genes.size() + "(" + filename + ") Genes with no hits:" + noHits );
        return resultRows;
    }

    public static Set<String> getUniqueGenes( Collection<String> rowNames ) {
        Set<String> genes = new HashSet<String>();
        for ( String row : rowNames ) {
            String geneName = row.replaceAll( "\\[.*\\]", "" );
            genes.add( geneName );
        }
        return genes;
    }

    public static double dotProduct( double[] a, double[] b ) {
        if ( a.length != b.length ) throw new RuntimeException( "Non matching array lengths" );
        double result = 0;
        for ( int i = 0; i < a.length; i++ ) {
            result += a[i] * b[i];
        }
        return result;
    }

    public static int intersectSize( Collection<?>... collections ) {
        return intersect( collections ).size();
    }

    public static Collection<?> subtract( Collection<?>... collections ) {
        Set<?> intersect = new HashSet<Object>( collections[0] );
        for ( int i = 1; i < collections.length; i++ ) {
            intersect.removeAll( collections[i] );
        }
        return intersect;
    }

    public static Collection<?> intersect( Collection<?>... collections ) {
        // make a new set so it doesnt modifiy the parameter sets
        Set<?> intersect = new HashSet<Object>( collections[0] );
        for ( int i = 1; i < collections.length; i++ ) {
            intersect.retainAll( collections[i] );
        }
        return intersect;
    }

    public static Collection<?> union( Collection<?>... collections ) {
        // make a new set so it doesnt modifiy the parameter sets
        Set<Object> union = new HashSet<Object>( collections[0] );
        for ( int i = 1; i < collections.length; i++ ) {
            union.addAll( collections[i] );
        }
        return union;
    }

    public static void addToVennMasterFile( File file, Collection items, String setName ) throws Exception {
        Collection<String> output = new LinkedList<String>();
        for ( Object item : items ) {
            output.add( "\"" + item.toString() + "\"\t" + setName );
        }
        FileTools.stringsToFile( output, file, true );
    }

    // public static void addToVennMasterFile( File file, Collection<String> items, String setName ) throws Exception {
    // Collection<String> output = new LinkedList<String>();
    // for ( Object item : items ) {
    // output.add( "\"" + item.toString() + "\"\t" + setName );
    // }
    // FileTools.stringsToFile( output, file, true );
    // }

    public static double zSum( DoubleMatrix<?, ?> x ) {
        DoubleMatrix2D a = DoubleFactory2D.dense.make( x.getRawMatrix() );
        return a.zSum();
    }
}
