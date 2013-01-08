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

import java.util.HashSet;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.Util;
import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;

public class IdentityAdjacency extends AdjacencyComputeSlow {
    private static Log log = LogFactory.getLog( IdentityAdjacency.class.getName() );

    public IdentityAdjacency( DoubleMatrix<String, String> matrix ) {
        super( matrix );
    }

    @Override
    public AdjacencyCompute clone() {
        // TODO Auto-generated method stub
        return this;
    }

    /**
     * Return the same matrix given.
     */
    public DoubleMatrix<String, String> getAdjacency( DoubleMatrix<String, String> toCompute ) {
        Set<String> rowSet = new HashSet<String>( toCompute.getRowNames() );
        Set<String> colSet = new HashSet<String>( toCompute.getColNames() );
        if ( !rowSet.equals( colSet ) ) {
            throw new RuntimeException( "Not a square matrix with matching rows/cols" );
        }
        // re-order the rows to be the same as cols
        DoubleMatrix<String, String> result = new DenseDoubleMatrix<String, String>( colSet.size(), colSet.size() );
        result.setColumnNames( toCompute.getColNames() );
        result.setRowNames( toCompute.getColNames() );
        for ( String row : colSet ) {
            for ( String col : colSet ) {
                result.setByKeys( row, col, toCompute.getByKeys( row, col ) );
                result.setByKeys( col, row, toCompute.getByKeys( row, col ) );
            }
        }
        return result;
    }

}
