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

import ubic.basecode.dataStructure.matrix.DoubleMatrix;

public interface AdjacencyCompute {
    public AdjacencyCompute clone();

    public void setupAdjacency();

    /**
     * The returned matrix must have the same column and row order as the input matrix column order
     * 
     * @param matrix
     * @param safe
     */
    public void setMatrix( DoubleMatrix<String, String> matrix, boolean safe );

    public DoubleMatrix<String, String> getAdjacency();

    public DoubleMatrix<String, String> getAdjacency( String removed, boolean adback );

    public DoubleMatrix<String, String> getAdjacency( String removed );

}
