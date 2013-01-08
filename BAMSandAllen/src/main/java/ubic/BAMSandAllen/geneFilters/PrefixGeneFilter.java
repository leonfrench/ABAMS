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
import java.util.Set;

import ubic.basecode.dataStructure.matrix.DoubleMatrix;

public class PrefixGeneFilter implements GeneFilter {
    List<String> prefixes;
    boolean removeMatches;

    public PrefixGeneFilter() {
        prefixes = new LinkedList<String>();
        removeMatches = true;
    }

    public String getName() {
        return "Prefix gene filter";
    }

    public PrefixGeneFilter( String prefix ) {
        this( prefix, false );
    }

    public PrefixGeneFilter( String prefix, boolean removeMatches ) {
        this();
        prefixes.add( prefix );
        this.removeMatches = removeMatches;
    }

    public List<String> getRowsToRemove( DoubleMatrix<String, String> matrix ) {
        List<String> rows = matrix.getRowNames();
        List<String> removeRows = new LinkedList<String>();

        for ( String rowName : rows ) {
            boolean matchesAPrefix = false;
            for ( String prefix : prefixes ) {
                if ( rowName.startsWith( prefix ) ) {
                    matchesAPrefix = true;
                }
            }
            // if we want matches removed and it matches, then remove it
            if ( matchesAPrefix && removeMatches ) removeRows.add( rowName );
            // if we want non matches removed, and it doesnt match, remove it
            if ( !matchesAPrefix && !removeMatches ) removeRows.add( rowName );
        }
        return removeRows;
    }
}
