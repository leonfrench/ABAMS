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

import ubic.BAMSandAllen.AllenDataLoaders.ImageSeriesInfoLoader;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;

public class PlaneRemoveFilter implements GeneFilter {
    public enum Plane {
        SAGITTAL, CORONAL
    };

    Plane plane;

    public PlaneRemoveFilter( Plane planeToRemove ) {
        this.plane = planeToRemove;
    }

    public String getName() {
        return "Plane remover, set to remove:" + plane + " image series";
    }

    public List<String> getRowsToRemove( DoubleMatrix<String, String> matrix ) {
        ImageSeriesInfoLoader loader = null;
        try {
            loader = new ImageSeriesInfoLoader();
        } catch ( Exception e ) {
            throw new RuntimeException( e );
        }
        List<String> rows = matrix.getRowNames();
        List<String> removeRows = new LinkedList<String>();

        for ( String rowName : rows ) {
            String planeString = loader.getPlaneFromRowName( rowName );
            if ( planeString.equals( "coronal" ) && plane.equals( Plane.CORONAL ) ) removeRows.add( rowName );
            if ( planeString.equals( "sagittal" ) && plane.equals( Plane.SAGITTAL ) ) removeRows.add( rowName );
        }
        return removeRows;
    }

}
