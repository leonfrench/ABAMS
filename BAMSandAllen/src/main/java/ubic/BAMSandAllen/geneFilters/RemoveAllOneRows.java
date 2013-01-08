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
package ubic.BAMSandAllen.geneFilters;

import java.util.LinkedList;
import java.util.List;

import ubic.BAMSandAllen.Util;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;

public class RemoveAllOneRows implements GeneFilter {

    
    //sorta redundant, could just use variance = 0 threshold
    public String getName() {
        return "All ones row remover";
    }

    public List<String> getRowsToRemove( DoubleMatrix<String, String> matrix ) {
        List<String> rows = matrix.getRowNames();
        List<String> removeRows = new LinkedList<String>();

        for ( String rowName : rows ) {
            double[] rowValues = matrix.getRowByName( rowName );
            boolean allOnes = false;
            for ( double value : rowValues ) {
                if ( value != 1 ) {
                    allOnes = false;
                    break;
                } else {
                    allOnes = true;
                }
            }
            if ( allOnes ) removeRows.add( rowName );
        }
        return removeRows;
    }
}
