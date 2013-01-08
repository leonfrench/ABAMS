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
import java.util.List;

import org.apache.commons.lang.time.StopWatch;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.RegressionVector;
import ubic.BAMSandAllen.SetupParameters;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenPartialExpressionMatrixPair;
import ubic.BAMSandAllen.MatrixPairs.MatrixPair;
import ubic.basecode.util.FileTools;

public class GreedyThreadRunner implements Runnable {
    private static Log log = LogFactory.getLog( GreedyThreadRunner.class.getName() );

    Collection<String> rowsToTest;
    MatrixPair pair;
    double baseLine;
    boolean increase;
    boolean slow;
    double bestIncrease;
    String bestRow;

    public GreedyThreadRunner( MatrixPair pair, boolean slow, boolean increase ) {
        this.slow = slow;
        this.pair = pair;
        this.increase = increase;
    }

    public void setRowsToTest( Collection<String> rowsToTest ) {
        this.rowsToTest = rowsToTest;
    }

    public void setBaseline( double baseLine ) {
        this.baseLine = baseLine;
    }

    public double getBestIncrease() {
        return bestIncrease;
    }

    public String getBestRow() {
        return bestRow;
    }

    public void removeRow( String row ) {
        pair.removeMatrixBDataRowFast( row );
    }

    // needs refactoring
    public void setTrianglesMatrixB( RegressionVector triangles ) {
        ( ( ConnectivityAndAllenPartialExpressionMatrixPair ) pair ).setTrianglesMatrixB( triangles );
    }

    public void run() {
        StopWatch watch = new StopWatch();
        watch.start();
        StopWatch smallWatch = new StopWatch();

        smallWatch.reset();
        smallWatch.start();
        // find the row that increases the correlation the most

        double baseCorrelation = baseLine;
        // log.info( "Base correlation:" + baseCorrelation + " size:" + pair.getMatrixBDataRows().size() );
        // pair.printDimensions();

        bestIncrease = Double.MAX_VALUE * -1;
        bestRow = null;
        int count = 0;

        if ( pair instanceof ConnectivityAndAllenPartialExpressionMatrixPair ) {
            ConnectivityAndAllenPartialExpressionMatrixPair partialPair = ( ( ConnectivityAndAllenPartialExpressionMatrixPair ) pair );
            // only do it if we are regressing on expression, ugly
            if ( !slow
                    && partialPair.getRegressType() == ConnectivityAndAllenPartialExpressionMatrixPair.RegressMatrix.BOTH
                    && partialPair.getRegressType() == ConnectivityAndAllenPartialExpressionMatrixPair.RegressMatrix.EXPRESSION ) {
                ( ( ConnectivityAndAllenPartialExpressionMatrixPair ) pair ).setRegressionComputeMatrixB( false );
            }
        }

        // if ( !slow && pair instanceof ConnectivityAndAllenPartialExpressionMatrixPair ) {
        // ( ( ConnectivityAndAllenPartialExpressionMatrixPair ) pair ).setRegressionComputeMatrixB( false );
        // }

        for ( String rowName : rowsToTest ) {
            // if ( ++count % 5000 == 0 ) log.info( "Data row:" + count + " time:" + watch.getTime() );

            double withoutCorrelation = pair.correWithoutMatrixBDataRow( rowName );
            double diff = withoutCorrelation - baseCorrelation;
            if ( !increase ) {
                diff *= -1;
            }

            // TEMP for ranking once
            // try {
            // FileTools.stringToFile( rowName + "\t" + diff + "\n", new File( SetupParameters.getDataFolder()
            // + "firstrun.txt" ), true );
            // } catch ( Exception e ) {
            // e.printStackTrace();
            // }
            // TEMP

            if ( diff > bestIncrease ) {
                bestRow = rowName;
                bestIncrease = diff;
            }
        }

        // if ( !slow && pair instanceof ConnectivityAndAllenPartialExpressionMatrixPair ) {
        // ( ( ConnectivityAndAllenPartialExpressionMatrixPair ) pair ).setRegressionComputeMatrixB( true );
        // }

        if ( pair instanceof ConnectivityAndAllenPartialExpressionMatrixPair ) {
            ConnectivityAndAllenPartialExpressionMatrixPair partialPair = ( ( ConnectivityAndAllenPartialExpressionMatrixPair ) pair );
            // only do it if we are regressing on expression, ugly
            if ( !slow
                    && ( partialPair.getRegressType() == ConnectivityAndAllenPartialExpressionMatrixPair.RegressMatrix.BOTH || partialPair
                            .getRegressType() == ConnectivityAndAllenPartialExpressionMatrixPair.RegressMatrix.EXPRESSION ) ) {
                ( ( ConnectivityAndAllenPartialExpressionMatrixPair ) pair ).setRegressionComputeMatrixB( false );
            }
        }

        // log.info( bestRow + " changes correlation by " + bestIncrease + " time:" + ( smallWatch.getTime() / 1000 )
        // + "s total:" + ( watch.getTime() / 1000 ) + "s" );

    }

    /**
     * @param args
     */
    public static void main( String[] args ) {
        // TODO Auto-generated method stub
    }

}
