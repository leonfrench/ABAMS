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
package ubic.BAMSandAllen.MatrixPairs;

import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.ABAMSDataMatrix;
import ubic.BAMSandAllen.AllenDataLoaders.AllenAtlasAnnotationLoader;
import ubic.BAMSandAllen.AllenDataLoaders.StructureCatalogLoader;
import ubic.BAMSandAllen.adjacency.AdjacencyCompute;
import ubic.BAMSandAllen.adjacency.EuclidAdjacency;
import ubic.BAMSandAllen.adjacency.LogEuclidAdjacency;

public class LiteratureSpaceMatrixPair extends LiteratureMatrixPair {
    private static Log log = LogFactory.getLog( LiteratureSpaceMatrixPair.class.getName() );

    public LiteratureSpaceMatrixPair( boolean logDistance ) throws Exception {
        // only ABA has space coordinates
        super( Nomenclature.ABA );

        try {
            AllenAtlasAnnotationLoader spaceLoader = new AllenAtlasAnnotationLoader();
            String matrixAname = "EuclidianDistance";

            AdjacencyCompute adjacencyCompute = null;
            if ( logDistance ) {
                adjacencyCompute = new LogEuclidAdjacency();
                matrixAname = "Log" + matrixAname;
            } else {
                adjacencyCompute = new EuclidAdjacency();
            }
            matrixA = new ABAMSDataMatrix( spaceLoader.getCenterMatrix(), matrixAname, adjacencyCompute );

            StructureCatalogLoader loader = new StructureCatalogLoader();
            Set<String> leafs = loader.getLeafs();
            leafs.retainAll( matrixA.getColNames() );
            // just the leafs
            matrixA = matrixA.retainColumns( leafs );

        } catch ( Exception e ) {
            e.printStackTrace();
            System.exit( 1 );
        }

    }

    /**
     * @param args
     */
    public static void main( String[] args ) throws Exception {
        boolean logDistance = false;
        LiteratureSpaceMatrixPair pair = new LiteratureSpaceMatrixPair( logDistance );
        pair.run();
        log.info( pair.getRegionsForWord( "macrosmatic" ));
        log.info( pair.getRegionsForWord( "488" ));
        log.info( pair.getRegionsForWord( "declining" ));
        log.info( pair.getRegionsForWord( "clarifies" ));
        log.info( pair.getRegionsForWord( "sexually" ));
        log.info( pair.getRegionsForWord( "from" ));
        log.info( pair.getRegionsForWord( "rostral" ));
        log.info( pair.getRegionsForWord( "projection" ));
    }

}
