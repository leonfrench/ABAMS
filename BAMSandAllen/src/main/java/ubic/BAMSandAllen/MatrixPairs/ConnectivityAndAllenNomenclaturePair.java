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

import java.util.HashSet;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.ABAMSDataMatrix;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes;
import ubic.BAMSandAllen.AllenDataLoaders.StructureCatalogLoader;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes.Direction;
import ubic.BAMSandAllen.ClassSelectors.BrainRegionClassSelector;
import ubic.BAMSandAllen.adjacency.DotProductAdjacency;

public class ConnectivityAndAllenNomenclaturePair extends ConnectivityAndAllenDataPair {
    private static Log log = LogFactory.getLog( ConnectivityAndAllenNomenclaturePair.class.getName() );
    boolean squareNomenclature;

    public ConnectivityAndAllenNomenclaturePair( BrainRegionClassSelector selector, boolean squareConnectivity,
            boolean squareNomenclature, Set<String> colNames, Direction direction ) throws Exception {
        super( selector, squareConnectivity, direction );
        this.squareNomenclature = squareNomenclature;
        // child or parent?

        StructureCatalogLoader loader = new StructureCatalogLoader();
        matrixB = new ABAMSDataMatrix( loader.getNomenclatureMatrix(), "Nomenclature", new DotProductAdjacency() );
        // log.info( matrixB );

        // reduce to the cols we were given
        if ( colNames != null ) {
            matrixB = matrixB.retainColumns( colNames );
        }
    }

    /*
     * can't do for parent or children when using leafs only. Might be usefull for sibling relationships
     */
    public void squareNomenclature() {
        log.info( "Squaring up" );
        matrixB = matrixB.getSquare();
        matrixA = removeZeroLines( matrixB );
    }

    public double run() throws Exception {
        super.run();
        if ( squareNomenclature ) {
            squareNomenclature();
        }
        return getCorrelation();
    }

    /**
     * @param args
     */
    public static void main( String[] args ) throws Exception {

        boolean squareMatrix = false;
        boolean virtualRegions = true;
        Direction direction = AnalyzeBAMSandAllenGenes.Direction.INCOMING;

        Set<String> colNames = new HashSet<String>( ExpressionMatrixPairFactory.getUsedCols( squareMatrix, direction ) );

        log.info( "________________________________________________________" );
        log.info( "________________________________________________________" );

        ConnectivityAndAllenNomenclaturePair x = new ConnectivityAndAllenNomenclaturePair(
                new BrainRegionClassSelector(), squareMatrix, false, colNames, direction );
        if ( virtualRegions ) x.makeVirtualRegions();
        x.run();
        x.runAllenStyle();

        // x.getCorrelation();
        // x.writeImages();
        // x.writeRMatrices();
        // x.test( 1000 );
    }
}
