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
package ubic.BAMSandAllen.AllenDataLoaders;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;

/**
 * Abstract class for a brain region loader, does not use ID's to represent regions, just strings
 * 
 * @author leon
 * 
 */
public abstract class RegionLoader {

    public int getRegionAmount() {
        return getRegions().size();
    }

    /*
     * get all regions, strings are exactly as seen in the source
     */
    public abstract Set<String> getRegions();

    /*
     * get all regions for lexicon purposes, trims and converts to lowercase
     */
    public Set<String> getRegionsForLexicon() {
        Set<String> unClean = getRegions();
        Set<String> result = new HashSet<String>();
        for ( String region : unClean ) {
            String newR = region.toLowerCase().trim();
            result.add( newR );
        }
        return result;
    }
    public DoubleMatrix<String, String> getNomenclatureMatrix() {
        boolean allParents = true;
        return getParentMatrix( allParents );
    }

    public DoubleMatrix<String, String> getParentMatrix() {
        boolean allParents = false;
        return getParentMatrix( allParents );
    }

    public DoubleMatrix<String, String> getParentMatrix( boolean allParents ) {

        List<String> regions = new LinkedList<String>( getRegions() );
        DoubleMatrix<String, String> result = new DenseDoubleMatrix<String, String>( regions.size(), regions.size() );
        result.setColumnNames( regions );
        result.setRowNames( regions );
        for ( String regionA : regions ) {
            System.out.println(regionA);
            for ( String regionB : regions ) {
                double value = 0d;
                // parents are columns, - entries in a regions column are all parents of it
                if ( ( allParents && getParents( regionA ).contains( regionB ) )
                        || ( !allParents && getParent( regionA ) != null && getParent( regionA ).contains( regionB ) ) ) {
                    // log.info( regionB + "->" + regionA);
                    value = 1d;
                }
                result.setByKeys( regionB, regionA, value );
            }
        }
        return result;
    }

    // get all leafs (string) - regions with no children
    public Set<String> getLeafs() {
        Set<String> leafs = new HashSet<String>();
        for ( String region : getRegions() ) {
            Set<String> children = getChildren( region );
            // if it has no children, then it's a leaf
            if ( children == null || children.size() == 0 ) {
                leafs.add( region );
            }
        }
        return leafs;
    }

    /*
     * Get how many leafs
     */
    public int getLeafAmount() {
        return getLeafs().size();
    }

    // get parent
    public abstract String getParent( String region );

    public abstract Set<String> getParents( String region );

    // get children (string), returns null if no children
    public abstract Set<String> getChildren( String region );
}
