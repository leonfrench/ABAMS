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
package ubic.BAMSandAllen.AllenDataLoaders.Explore;

import ubic.BAMSandAllen.AllenDataLoaders.StructureCatalogLoader;

public class ShowParents {
    StructureCatalogLoader catalog;

    public ShowParents() throws Exception {
        catalog = new StructureCatalogLoader();
    }

    public void printParents( String region ) {
        if ( region == null ) return;
        printParents( catalog.getParent( region ) );
        System.out.print( "," + region );
    }

    public int countParents( String region ) {
        if ( region == null ) return -1;
        return 1 + countParents( catalog.getParent( region ) );
    }

    public void printAllParents() {
        for ( String region : catalog.getRegions() ) {
            System.out.print( region );
            System.out.print(","+countParents( region ));
            printParents( region );
            System.out.println();
        }
    }

    /**
     * @param args
     */
    public static void main( String[] args ) throws Exception {
        // TODO Auto-generated method stub
        // get all the regions and show the heirachy
        ShowParents showParents = new ShowParents();
        showParents.printAllParents();
    }
}
