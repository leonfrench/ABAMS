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
package ubic.BAMSandAllen.adjacency;

import javax.vecmath.Point3d;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenExpressionMatrixPair;
import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;

public class ChiIndexAdjacency extends AdjacencyComputeSlow implements AdjacencyCompute {
    private static Log log = LogFactory.getLog( ChiIndexAdjacency.class.getName() );

    // from http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1676027/bin/pcbi.0020167.sd001.pdf
    // Gene Expression of Caenorhabditis elegans Neurons Carries Information on Their Synaptic Connectivity
    // Alon Kaufman, Gideon Dror, Isaac Meilijson, and Eytan Ruppin

    public DoubleMatrix<String, String> getAdjacency( DoubleMatrix<String, String> input ) {
        int size = input.columns();
        DoubleMatrix<String, String> result = new DenseDoubleMatrix<String, String>( size, size );
        result.setColumnNames( input.getColNames() );
        result.setRowNames( input.getColNames() );

        for ( int i = 0; i < size; i++ ) {
            double[] iCol = input.getColumn( i );
            for ( int j = i; j < size; j++ ) {
                double a, b, c, d;
                a = b = c = d = 0;
                double[] jCol = input.getColumn( j );
                for ( int rowIndex = 0; rowIndex < jCol.length; rowIndex++ ) {
                    boolean x = 1d == iCol[rowIndex];
                    boolean y = 1d == jCol[rowIndex];
                    // if ( Double.isNaN( x ) || Double.isNaN( y ) ) continue;
                    if ( !x && !y ) a++;
                    if ( x && !y ) b++;
                    if ( !x && y ) c++;
                    if ( x && y ) d++;
                }
                double n = a + b + c + d;
                // log.info( a );
                // log.info( b );
                // log.info( c );
                // log.info( d );
                if ( n != jCol.length ) throw new RuntimeException( "Error in chisquared, sums dont match, NaNs?" );
                double chiSq = ( n * Math.pow( a * d - c * b, 2 ) ) / ( ( a + b ) * ( c + d ) * ( a + c ) * ( b + d ) );
                double chi = Math.sqrt( n / ( ( a + b ) * ( c + d ) * ( a + c ) * ( b + d ) ) )
                        * Math.abs( a * d - c * b );
                // chi2 = math.pow(chi,2)
                result.set( i, j, chiSq );
                result.set( j, i, chiSq );
            }
        }
        return result;
    }

    public ChiIndexAdjacency clone() {
        return this;
    }
}
