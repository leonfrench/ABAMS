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

import java.io.Serializable;
import java.util.LinkedList;
import java.util.List;

import ubic.BAMSandAllen.Util;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;

public abstract class AdjacencyComputeSlow implements AdjacencyCompute, Serializable {
    DoubleMatrix<String, String> matrix;
    DoubleMatrix<String, String> cachedAdjacency;

    public AdjacencyComputeSlow( DoubleMatrix<String, String> matrix ) {
        setMatrix( matrix, true );
    }

    public AdjacencyComputeSlow() {
        // expect a setmatrix call
        this.matrix = null;
    }

    public final DoubleMatrix<String, String> getAdjacency() {
        if ( cachedAdjacency == null ) {
            cachedAdjacency = getAdjacency( matrix );
        }
        return cachedAdjacency;
    }

    public abstract DoubleMatrix<String, String> getAdjacency( DoubleMatrix<String, String> toCompute );

    public DoubleMatrix<String, String> getAdjacency( String removed, boolean adback ) {
        throw new RuntimeException( "not supported" );
    }

    public DoubleMatrix<String, String> getAdjacency( String removed ) {
        List<String> rows = new LinkedList<String>();
        DoubleMatrix<String, String> matrixCompute = Util.removeRows( matrix, rows );
        return getAdjacency( matrixCompute );
    }

    public void setMatrix( DoubleMatrix<String, String> matrix, boolean safe ) {
        this.matrix = matrix;
        cachedAdjacency = null;
    }

    public abstract AdjacencyCompute clone();

    public void setupAdjacency() {
        // no setup
    }

}
