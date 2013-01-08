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

import java.io.FileInputStream;
import java.io.ObjectInputStream;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.ABAMSDataMatrix;
import ubic.BAMSandAllen.SetupParameters;
import ubic.BAMSandAllen.StructureCatalogAnalyze;
import ubic.BAMSandAllen.Util;
import ubic.BAMSandAllen.AllenDataLoaders.StructureCatalogLoader;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes.Direction;
import ubic.BAMSandAllen.ClassSelectors.BrainRegionClassSelector;
import ubic.BAMSandAllen.adjacency.CorrelationAdjacency;
import ubic.BAMSandAllen.geneFilters.RemoveAllOneRows;
import ubic.BAMSandAllen.geneFilters.RemoveLowVarianceGeneFilter;
import ubic.BAMSandAllen.optimize.GreedyMultiThreaded;
import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;

public class ConnectivityLiteratureMatrixPair extends LiteratureMatrixPair {
    private static Log log = LogFactory.getLog( ConnectivityLiteratureMatrixPair.class.getName() );

    boolean squareConnectivity;
    Direction direction;

    public ConnectivityLiteratureMatrixPair( Direction direction ) throws Exception {
        super( Nomenclature.BAMS );

        // allenCatalog = new StructureCatalogLoader();
        squareConnectivity = false;
        this.direction = direction;
        BrainRegionClassSelector selector = new BrainRegionClassSelector();

        if ( squareConnectivity && direction.equals( Direction.APPENDED ) ) {
            throw new RuntimeException( "Error square matrix incompatable with appended connectivity" );
        }
        StructureCatalogAnalyze forMatrix = new StructureCatalogAnalyze( selector );
        // could be an option here for propigated or non-propigated connection matrix
        forMatrix.readModel( SetupParameters.getDataFolder() + "Propigated.rdf" );

        DoubleMatrix<String, String> dataMatrix = forMatrix.makeConnectionMatrix( direction );

        matrixA = new ABAMSDataMatrix( dataMatrix, "Connectivity", new CorrelationAdjacency( dataMatrix ) );

        // get Allen leaf regions and convert them to BAMS regions
        StructureCatalogLoader loader = new StructureCatalogLoader();
        Collection<String> leafsABA = loader.getLeafs();
        Collection<String> leafsBAMS = new HashSet<String>();
        for ( String leafABA : leafsABA ) {
            leafsBAMS.addAll( loader.getBAMSMappedRegions( leafABA ) );
        }

        // remove the non leaf nodes
        Collection<String> nonLeafs = ( Collection<String> ) Util.subtract( matrixA.getColNames(), leafsBAMS );
        log.info( "Removing " + nonLeafs.size() + " regions " + nonLeafs.toString() );
        matrixA = matrixA.removeColumns( nonLeafs );

        // log.info( "Removing bed Nuclei Stria" );
        // removeBedNucleiStria();

        // matrixA = matrixA.removeRows( Util.findZeroRows( matrixA ) );

        matrixA = matrixA.removeZeroColumns();
        matrixB = matrixB.removeZeroColumns();

        log.info( "got Matrix A - connections" );
    }


    /**
     * @param args
     */
    public static void main( String[] args ) throws Exception {
        ConnectivityLiteratureMatrixPair pair = new ConnectivityLiteratureMatrixPair( Direction.INCOMING );
        pair.printDimensions();
        pair.run();
        log.info( pair.getMatrixB().getRowByName( "the" ) );
        Set<String> rows = new HashSet<String>( pair.getMatrixB().getRowNames() );
        pair.removeLiteratureOnes();
        pair.printDimensions();

        log.info( "subtracted:" + Util.subtract( rows, pair.getMatrixB().getRowNames() ) );
        // pair.writeImages();

        // System.exit( 1 );
        int threads = 16;
        GreedyMultiThreaded remover = new GreedyMultiThreaded( pair, threads );
        boolean keepSign = true; // should it increase or decrease the correlation?
        remover.run( 22076, keepSign );
    }

}
