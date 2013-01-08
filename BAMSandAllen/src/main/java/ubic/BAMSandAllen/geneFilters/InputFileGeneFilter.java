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

public class InputFileGeneFilter implements GeneFilter {
    String filename;

    public InputFileGeneFilter( String filename ) {
        this.filename = filename;
    }

    public List<String> getRowsToRemove( DoubleMatrix<String, String> matrix ) {
        try {
            List<String> fromFile = Util.getRowNamesFromGeneNames( filename, matrix.getRowNames() );
            List<String> result = new LinkedList<String>( matrix.getRowNames() );
            result.removeAll( fromFile );
            return result;
        } catch ( Exception e ) {
            throw new RuntimeException( e );
        }

    }

    public String getName() {
        return "Input file filter from:" + filename;
    }
}
