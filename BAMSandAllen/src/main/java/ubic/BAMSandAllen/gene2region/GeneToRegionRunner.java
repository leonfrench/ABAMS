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
package ubic.BAMSandAllen.gene2region;

import java.io.File;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.Util;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes.Direction;
import ubic.BAMSandAllen.MatrixPairs.AllenMatrixPairFactory;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenExpressionMatrixPair;
import ubic.BAMSandAllen.MatrixPairs.ExpressionMatrixPairFactory;
import ubic.BAMSandAllen.MatrixPairs.MatrixPair;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenPartialExpressionMatrixPair.RegressMatrix;
import ubic.BAMSandAllen.geneFilters.NaNGeneFilter;
import ubic.BAMSandAllen.geneFilters.NonExpFilter;
import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;

public class GeneToRegionRunner {
    private static Log log = LogFactory.getLog( GeneToRegionRunner.class.getName() );

    ConnectivityAndAllenExpressionMatrixPair pair;

    public GeneToRegionRunner( ConnectivityAndAllenExpressionMatrixPair pair ) {
        this.pair = pair;
    }

    public boolean hasNonZeroEntry( double[] ar ) {
        for ( double d : ar ) {
            if ( d != 0 ) return true;
        }
        return false;
    }

    public Set<String> getNonNaNRows( DoubleMatrix<String, String> matrix ) {
        Set<String> result = new HashSet<String>();
        for ( String row : matrix.getRowNames() ) {
            for ( double d : matrix.getRowByName( row ) ) {
                if ( Double.isNaN( d ) ) {
                    result.add( row );
                    break;
                }
            }
        }
        return result;
    }

    public void runPartialAll( DoubleMatrix<String, String> explain, boolean jacknife, RegressMatrix regressType )
            throws Exception {
        DoubleMatrix<String, String> expression = pair.getMatrixB();
        DoubleMatrix<String, String> connectivity = pair.getMatrixA();

        List<String> nonNaNRows = new LinkedList<String>( getNonNaNRows( expression ) );
        log.info( "Non nan rows:" + nonNaNRows.size() + " of " + expression.rows() );

        DoubleMatrix<String, String> result = new DenseDoubleMatrix<String, String>( expression.rows(), expression
                .columns() );
        // limited to ABA mapped regions
        result.setRowNames( expression.getRowNames() );
        result.setColumnNames( expression.getColNames() );

        File f = new File( pair.getBaseName() + ".GeneToRegion.txt" );
        int count = 0;
        GeneToRegionCreator gene2region;

        // use limited set of regions
        for ( String region : expression.getColNames() ) {
            if ( !connectivity.containsRowName( region ) ) continue;
            for ( String gene : expression.getRowNames() ) {
                gene2region = new GeneToRegionCreator( gene, region, expression, connectivity, explain, regressType );
                if ( jacknife ) {
                    result.setByKeys( gene, region, gene2region.getJackKnifeCorrelation() );
                } else {
                    result.setByKeys( gene, region, gene2region.getCorrelation() );
                }
                // FileTools.stringToFile( gene2region.getCorrelation() + "\n", f, true );
                // FileTools.stringToFile( gene2region.toString() + "\n", f, true );

                if ( ++count % 10000 == 0 ) log.info( count + " " + gene2region.getCorrelation() );
                // log.info( gene2region.getCorrelation() );
            }
            // System.exit( 1 );
        }
        if ( jacknife ) {
            Util.writeRTable( pair.getBaseName() + ".GeneToRegion.partial.exp.jack.txt", result );
        } else {
            Util.writeRTable( pair.getBaseName() + ".GeneToRegion.partial.exp.txt", result );
        }
    }

    public void runAll() throws Exception {
        DoubleMatrix<String, String> expression = pair.getMatrixB();
        DoubleMatrix<String, String> connectivity = pair.getMatrixA();

        DoubleMatrix<String, String> result = new DenseDoubleMatrix<String, String>( expression.rows(), connectivity
                .rows() );
        result.setRowNames( expression.getRowNames() );
        result.setColumnNames( connectivity.getRowNames() );

        // File f = new File( pair.getBaseName() + ".GeneToRegion.txt" );
        int count = 0;
        GeneToRegionCreator gene2region;

        for ( String region : connectivity.getRowNames() ) {

            for ( String gene : expression.getRowNames() ) {
                gene2region = new GeneToRegionCreator( gene, region, expression, connectivity );
                result.setByKeys( gene, region, gene2region.getJackKnifeCorrelation() );
                // FileTools.stringToFile( gene2region.getCorrelation() + "\n", f, true );
                // FileTools.stringToFile( gene2region.toString() + "\n", f, true );

                if ( ++count % 10000 == 0 ) log.info( count );
                // log.info( gene2region.getCorrelation() );
            }
        }
        Util.writeRTable( pair.getBaseName() + ".GeneToRegion.jacknife.txt", result );
    }

    /**
     * @param args
     */
    public static void main( String[] args ) throws Exception {
        // TODO Auto-generated method stub
        boolean useVirtual = true;
        boolean removeNonExp = true;
        boolean logDistance = true;
        Direction direction = Direction.APPENDED;
        ConnectivityAndAllenExpressionMatrixPair full = ExpressionMatrixPairFactory.connectivityAndExpression(
                direction, useVirtual, removeNonExp );
        full.applyGeneFilter( new NaNGeneFilter( 30 ) );
        full.printDimensions();
        log.info( full.getCorrelation() );
        // System.exit( 1 );

        // runner.runAll();

        log.info( "______________________________" );
        log.info( "______________________________" );
        GeneToRegionRunner runner = new GeneToRegionRunner( full );

        runner.runAll();
        System.exit( 1 );

        MatrixPair explainPair = AllenMatrixPairFactory.getConnectivityAndDistancePair( direction, useVirtual,
                logDistance );
        runner.runPartialAll( explainPair.getAdjacencyB(), false, RegressMatrix.EXPRESSION );
        runner.runPartialAll( explainPair.getAdjacencyB(), true, RegressMatrix.EXPRESSION );
    }
}
