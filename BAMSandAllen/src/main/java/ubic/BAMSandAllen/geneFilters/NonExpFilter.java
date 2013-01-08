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

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.SetupParameters;
import ubic.BAMSandAllen.TopTenInfo;
import ubic.BAMSandAllen.Util;
import ubic.BAMSandAllen.AllenDataLoaders.ImageSeriesInfoLoader;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;

public class NonExpFilter implements GeneFilter {
    private static Log log = LogFactory.getLog( TopTenInfo.class.getName() );
    boolean removeCoronal;

    public NonExpFilter() {
        this( false );
    }

    public NonExpFilter( boolean removeCoronal ) {
        this.removeCoronal = removeCoronal;
    }

    public String getName() {
        // TODO Auto-generated method stub
        return "Allen Brain atlas Non expressed genes list";
    }

    public List<String> getRowsToRemove( DoubleMatrix<String, String> matrix ) {
        String filename = SetupParameters.getDataFolder() + "ABANonexpressed.txt";
        try {
            ImageSeriesInfoLoader imageInfo = new ImageSeriesInfoLoader();
            List<String> nonExp = Util.getRowNamesFromGeneNames( filename, matrix.getRowNames() );
            List<String> result = new LinkedList<String>();
            int counter = 0;
            for ( String row : nonExp ) {
                // if it does not have a coronal section, add it - resolves disagreement between AGEA gene set and non
                // expressing list
                if ( !removeCoronal || !imageInfo.hasCoronalImageFromRowName( row ) )
                    result.add( row );
                else {
                    counter++;
                }
            }
            log.info( counter + " image series sets were removed from non exp list for having a cornal imageset" );
            return result;
        } catch ( Exception e ) {
            throw new RuntimeException( e );
        }
    }
}
