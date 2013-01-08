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

import ubic.BAMSandAllen.SetupParameters;
import ubic.BAMSandAllen.Util;
import ubic.BAMSandAllen.AllenDataLoaders.ImageSeriesInfoLoader;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.util.FileTools;

public class AGEAGeneFilter implements GeneFilter {
    String filename;
    boolean imageSeries;

    public AGEAGeneFilter( boolean imageSeries ) {
        this.imageSeries = imageSeries;
        if ( imageSeries ) {
            filename = SetupParameters.getDataFolder() + "AGEA Imagesets.txt";
        } else {
            filename = SetupParameters.getDataFolder() + "AGEA genes.txt";
        }
    }

    public String getName() {
        // TODO Auto-generated method stub
        return "Allen Brain atlas AGEA gene list. Filename:" + filename;
    }

    public List<String> getRowsToRemove( DoubleMatrix<String, String> matrix ) {

        if ( imageSeries == false ) {
            try {
                List<String> keepRows = Util.getRowNamesFromGeneNames( filename, matrix.getRowNames() );
                List<String> result = new LinkedList<String>( matrix.getRowNames() );
                result.removeAll( keepRows );
                return result;
            } catch ( Exception e ) {
                throw new RuntimeException( e );
            }
        } else {
            List<String> result = new LinkedList<String>();
            try {
                List<String> imageSeriesSets = FileTools.getLines( filename );
                for ( String row : matrix.getRowNames() ) {
                    String imageID = ImageSeriesInfoLoader.getImageIDFromRowLabel( row );
                    if ( imageSeriesSets.contains( imageID ) ) {
                        // then keep it
                    } else {
                        result.add( row );
                    }
                }
                return result;
            } catch ( Exception e ) {
                throw new RuntimeException( e );
            }
        }
    }
}
