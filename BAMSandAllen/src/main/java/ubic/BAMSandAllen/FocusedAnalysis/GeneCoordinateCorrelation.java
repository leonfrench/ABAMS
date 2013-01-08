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

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.ABAMSDataMatrix;
import ubic.BAMSandAllen.Util;
import ubic.BAMSandAllen.ClassSelectors.BrainRegionClassSelector;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenExpressionMatrixPair;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenSpacePair;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;


/**
 * Compares gene expression levels to various brain region statistics like position and connectivity degree.
 * 
 * @author leon
 *
 */
public class GeneCoordinateCorrelation {
    private static Log log = LogFactory.getLog( ConnectivityAndAllenExpressionMatrixPair.class.getName() );

    ConnectivityAndAllenExpressionMatrixPair expressionPair;
    DoubleMatrix<String, String> spaceMatrix;

    public GeneCoordinateCorrelation( ConnectivityAndAllenExpressionMatrixPair pair ) throws Exception {
        expressionPair = pair;
        // get a BAMS version of the location matrix
        ConnectivityAndAllenSpacePair spacePair = new ConnectivityAndAllenSpacePair( new BrainRegionClassSelector(),
                pair.isSquareConnectivity(), null, pair.getDirection() );

        if ( pair.hasVirtualRegions() ) spacePair.makeVirtualRegions();
        spacePair.run();
        spaceMatrix = spacePair.getMatrixB();

    }

    // done in BAMS space
    public double getDegreeCorrelation( String axis ) {
        ABAMSDataMatrix matrixA = expressionPair.getMatrixA();
        double[] degrees = Util.columnSums( matrixA ).getRow( 0 );
        double[] coordinates = new double[degrees.length];
        for ( int i = 0; i < coordinates.length; i++ ) {
            coordinates[i] = spaceMatrix.getByKeys( axis, matrixA.getColName( i ) );
        }
        return Util.spearmanCorrel( coordinates, degrees );
    }

    public double getX( String colname ) {
        return spaceMatrix.getByKeys( "x", colname );
    }

    public double getY( String colname ) {
        return spaceMatrix.getByKeys( "y", colname );
    }

    public double getZ( String colname ) {
        return spaceMatrix.getByKeys( "z", colname );
    }

    public double getXRankCorrelation( String rowName ) {
        ABAMSDataMatrix matrixB = expressionPair.getMatrixB();
        double[] exp = matrixB.getRowByName( rowName );
        double[] coordinates = new double[exp.length];
        for ( int i = 0; i < coordinates.length; i++ ) {
            coordinates[i] = spaceMatrix.getByKeys( "x", matrixB.getColName( i ) );
        }
        return Util.spearmanCorrel( coordinates, exp );
    }

    public double getYRankCorrelation( String rowName ) {
        ABAMSDataMatrix matrixB = expressionPair.getMatrixB();
        double[] exp = matrixB.getRowByName( rowName );
        double[] coordinates = new double[exp.length];

        for ( int i = 0; i < coordinates.length; i++ ) {
            coordinates[i] = spaceMatrix.getByKeys( "y", matrixB.getColName( i ) );
        }
        return Util.spearmanCorrel( coordinates, exp );
    }

    public double getZRankCorrelation( String rowName ) {
        ABAMSDataMatrix matrixB = expressionPair.getMatrixB();
        double[] exp = matrixB.getRowByName( rowName );
        double[] coordinates = new double[exp.length];

        for ( int i = 0; i < coordinates.length; i++ ) {
            coordinates[i] = spaceMatrix.getByKeys( "z", matrixB.getColName( i ) );
        }
        return Util.spearmanCorrel( coordinates, exp );
    }

}
