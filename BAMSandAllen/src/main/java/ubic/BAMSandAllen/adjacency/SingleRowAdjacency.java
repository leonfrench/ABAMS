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

import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;

public class SingleRowAdjacency extends AdjacencyComputeSlow {
    String rowName;
    boolean add;

    public SingleRowAdjacency( String rowName, boolean add ) {
        super();
        this.rowName = rowName;
        this.add = add;
    }

    /**
     * @param matrix - should have only one row
     */
    public SingleRowAdjacency( DoubleMatrix<String, String> matrix, String rowName, boolean add ) {
        super( matrix );
        this.rowName = rowName;
        this.add = add;
    }

    public AdjacencyCompute clone() {
        return this;
    }

    public DoubleMatrix<String, String> getAdjacency( DoubleMatrix<String, String> input ) {
        int size = input.columns();
        DoubleMatrix<String, String> result = new DenseDoubleMatrix<String, String>( size, size );
        result.setColumnNames( input.getColNames() );
        result.setRowNames( input.getColNames() );
        int rowIndex = input.getRowIndexByName( rowName );
        for ( int i = 0; i < size; i++ ) {
            double iValue = input.get( rowIndex, i );
            for ( int j = i; j < size; j++ ) {
                double jValue = input.get( rowIndex, j );
                // double value = Math.log( iValue / jValue );
                double value;
                // awkward
                if ( !add )
                    value = iValue - jValue;
                else
                    value = iValue + jValue;

                result.set( i, j, value );
                result.set( j, i, value );
            }
        }
        return result;
    }

}
