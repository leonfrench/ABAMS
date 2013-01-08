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
package ubic.BAMSandAllen.FocusedAnalysis;

import java.io.FileWriter;
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.commons.math.stat.StatUtils;
import org.apache.commons.math.stat.inference.TestUtils;

import com.sun.org.apache.xalan.internal.xsltc.compiler.sym;

import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes;
import ubic.BAMSandAllen.SetupParameters;
import ubic.BAMSandAllen.Util;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes.Direction;
import ubic.BAMSandAllen.BAMSDataLoaders.BAMSDataLoader;
import ubic.BAMSandAllen.MappingHelpers.BamsAllenMapper;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenExpressionMatrixPair;
import ubic.BAMSandAllen.MatrixPairs.ExpressionMatrixPairFactory;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenPartialExpressionMatrixPair.RegressMatrix;
import ubic.BAMSandAllen.adjacency.IdentityAdjacency;
import ubic.BAMSandAllen.geneFilters.NaNGeneFilter;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;

public class SquareTests {
    private static Log log = LogFactory.getLog( SquareTests.class.getName() );

    /**
     * @param args
     */
    public static void main( String[] args ) throws Exception {
        Direction direction = AnalyzeBAMSandAllenGenes.Direction.ANYDIRECTION;

        ConnectivityAndAllenExpressionMatrixPair pair;
        

        boolean logDistance = true;
        boolean slow = false;
        RegressMatrix regressType = RegressMatrix.BOTH;
        // cannot use virtual regions because it will lead to a non-square connectivity matrix (incoming only) - virtual
        // regions only formed on the columns, not rows - fixed
        boolean useVirtual = true;
        boolean removeNonExp = true;
        boolean run = true;
        boolean squareMatrix = true;

        direction = AnalyzeBAMSandAllenGenes.Direction.ANYDIRECTION;
        pair = ExpressionMatrixPairFactory.connectivityAndExpression( direction, useVirtual, removeNonExp, run,
                squareMatrix );

        pair = make142SquarePair( direction );
        pair.writeRMatrices();
        
        

        //log.info( "Result:" + testOldHypoth( pair ) );
        //System.exit( 1 );

        direction = AnalyzeBAMSandAllenGenes.Direction.INCOMING;

        pair = ExpressionMatrixPairFactory.connectivityAndExpression( direction, useVirtual, removeNonExp, run,
                squareMatrix );
        pair.printDimensions();
        // pair.squareConnectionMatrix();
        pair.getMatrixA().setAdjacencyCompute( new IdentityAdjacency( pair.getMatrixA() ) );
        // pair.writeRMatrices();
        // pair.writeImages();
        System.out.println( pair.getCorrelation() );
        //System.exit( 1 );

        pair.test( 1000, false );

        pair = ExpressionMatrixPairFactory.connectivityPartial( direction, slow, regressType, useVirtual, removeNonExp,
                logDistance, squareMatrix );
        pair.test( 1000, false );

    }

    public static ConnectivityAndAllenExpressionMatrixPair make142SquarePair( Direction direction ) throws Exception {
        boolean useVirtual = true;
        boolean removeNonExp = true;
        boolean run = false;
        boolean squareMatrix = true;

        ConnectivityAndAllenExpressionMatrixPair pair;

        pair = ExpressionMatrixPairFactory.connectivityAndExpression( direction, useVirtual, removeNonExp, run,
                squareMatrix );
        pair.printDimensions();

        boolean removeZeroLines = false;

        pair.slimMatrices();
        pair.testZeroes();
        pair.printDimensions();
        pair.removeConnectionZeroes();
        // if ( squareMatrix ) pair.squareConnectionMatrix();
        // if ( squareMatrix ) pair.squareConnectionMatrix();
        // pair.printMappingRelations();
        //
        // pair.printDimensions();
        // pair.testZeroes();
        //
        pair.squareConnectionMatrix( removeZeroLines );
        pair.sameSpace();
        // //pair.run();
        
       
        // for some reason it doesnt remove the NaN genes at the right point
        pair.applyGeneFilter( new NaNGeneFilter() );
        pair.printDimensions();

        return pair;

    }

    public static double testOldHypoth( ConnectivityAndAllenExpressionMatrixPair pair ) throws Exception {

        boolean symetric = pair.getDirection().equals( Direction.ANYDIRECTION );
        log.info( "Symetric:" + symetric );
        // matrix B is expression
        // DoubleMatrix<String, String> bCorrelations = Util.correlateColumns( matrixB, spearman );

        DoubleMatrix<String, String> bCorrelations = pair.getMatrixB().getAdjacency();
        DoubleMatrix<String, String> connectionMatrix = pair.getMatrixA();

        if ( bCorrelations.getRowNames().size() != connectionMatrix.getRowNames().size() ) {
            log.info( "Error, row size does not match" );
            return -1;
        }
        if ( bCorrelations.getColNames().size() != connectionMatrix.getColNames().size() ) {
            log.info( "Error, column size does not match" );
            return -1;
        }

        List<Double> connectedCorrelations = new LinkedList<Double>();
        List<Double> disconnectedCorrelations = new LinkedList<Double>();

        List<String> regionList = bCorrelations.getColNames();
        for ( int i = 0; i < regionList.size(); i++ ) {
            String regionOne = regionList.get( i );
            for ( int j = 0; j < regionList.size(); j++ ) {
                // only do one half of the matrix if it is symetric
                if ( symetric && i > j ) continue;
                String regionTwo = regionList.get( j );
                // skip diagonal
                if ( regionOne.equals( regionTwo ) ) continue;

                double correlation = bCorrelations.getByKeys( regionOne, regionTwo );
                double connected = connectionMatrix.getByKeys( regionOne, regionTwo );
                if ( connected == 1d ) {
                    connectedCorrelations.add( correlation );
                } else {
                    disconnectedCorrelations.add( correlation );
                }

            }
        }

        // log.info( disconnectedCorrelations );
        // log.info( connectedCorrelations );
        log.info( "Disconnected size:" + disconnectedCorrelations.size() );

        log.info( "Connected size:" + connectedCorrelations.size() );
        double[] discon = convert( disconnectedCorrelations );
        double[] conn = convert( connectedCorrelations );

        log.info( "Disconnected mean:" + StatUtils.mean( discon ) );
        log.info( "Variance:" + StatUtils.variance( discon ) );
        log.info( "Connected mean:" + StatUtils.mean( conn ) );
        log.info( "Variance:" + StatUtils.variance( conn ) );

        FileWriter f = new FileWriter( SetupParameters.getDataFolder() + "disconnect.txt" );
        f.write( disconnectedCorrelations.toString().replace( ']', ' ' ).replace( '[', ' ' ) );
        f.close();

        f = new FileWriter( SetupParameters.getDataFolder() + "connect.txt" );
        f.write( connectedCorrelations.toString().replace( ']', ' ' ).replace( '[', ' ' ) );
        f.close();

        return TestUtils.tTest( discon, conn );
    }

    public static double[] convert( List<Double> in ) {
        double[] out = new double[in.size()];
        for ( int i = 0; i < out.length; i++ )
            out[i] = in.get( i );
        return out;
    }

}
