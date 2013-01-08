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

import ubic.BAMSandAllen.Util;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;

public class LogEuclidAdjacency extends EuclidAdjacency {

    public LogEuclidAdjacency() {
        super();
    }

    public LogEuclidAdjacency( DoubleMatrix<String, String> matrix ) {
        super( matrix );
    }

    public LogEuclidAdjacency clone() {
        return this;
    }

    public DoubleMatrix<String, String> getAdjacency( DoubleMatrix<String, String> input ) {
        DoubleMatrix<String, String> euclidAdjacency = super.getAdjacency( input );
        // only in the case of one to many mappings will zero occur. When using virtual regions there should be no 0
        // distances
        DoubleMatrix<String, String> result = Util.logMatrix( euclidAdjacency, 1 );
        return result;
    }
}
