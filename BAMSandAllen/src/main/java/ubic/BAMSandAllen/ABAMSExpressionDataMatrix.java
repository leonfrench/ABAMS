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

import java.util.LinkedList;
import java.util.List;

import ubic.BAMSandAllen.adjacency.AdjacencyCompute;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;

public class ABAMSExpressionDataMatrix extends ABAMSDataMatrix {

    public ABAMSExpressionDataMatrix( DoubleMatrix<String, String> matrix, String name, AdjacencyCompute adjacency ) {
        super( matrix, name, adjacency );
    }

    public ABAMSExpressionDataMatrix( int rows, int cols, String name, AdjacencyCompute adjacency ) {
        super( rows, cols, name, adjacency );
    }

    public List<String> getNullNameRows() {
        List<String> removeRows = new LinkedList<String>();
        for ( String rowName : getRowNames() ) {
            if ( rowName.startsWith( "null[" ) || rowName.startsWith( "none[" ) ) removeRows.add( rowName );
        }
        return removeRows;
    }

}
