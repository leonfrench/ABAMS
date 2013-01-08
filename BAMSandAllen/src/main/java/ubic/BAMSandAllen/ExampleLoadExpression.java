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
package ubic.BAMSandAllen;

import java.util.Collection;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.AllenDataLoaders.AllenCatalogMatrices2;
import ubic.BAMSandAllen.AllenDataLoaders.StructureCatalogLoader;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenDataPair;
import ubic.BAMSandAllen.MatrixPairs.ExpressionMatrixPairFactory;
import ubic.BAMSandAllen.adjacency.CorrelationAdjacency;
import ubic.BAMSandAllen.geneFilters.GeneFilter;
import ubic.BAMSandAllen.geneFilters.NaNGeneFilter;
import ubic.BAMSandAllen.geneFilters.PlaneRemoveFilter;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;

public class ExampleLoadExpression {
    private static Log log = LogFactory.getLog( ExampleLoadExpression.class.getName() );

    public static void getAll() throws Exception {
        AllenCatalogMatrices2 allenMatrices = new AllenCatalogMatrices2();
        DoubleMatrix<String, String> dataMatrix = allenMatrices.getFromDisk( "NewEnergies" );
        ABAMSDataMatrix matrixB = new ABAMSDataMatrix( dataMatrix, "Energy", new CorrelationAdjacency( dataMatrix ) );
        Util.writeRTable( SetupParameters.getDataFolder() + "FullMatrix.txt", matrixB );

    }

    public static void getAllenProcessed() throws Exception {
        AnalyzeBAMSandAllenGenes.Direction direction = AnalyzeBAMSandAllenGenes.Direction.ANYDIRECTION;
        boolean removeNonExp = true;
        ABAMSDataMatrix geneRegionExpressionMatrix = ExpressionMatrixPairFactory.getEnergyMatrix( direction,
                removeNonExp );
        Util.writeRTable( SetupParameters.getDataFolder() + "ProcessedMatrix.txt", geneRegionExpressionMatrix );
    }

    public static void getABAMSProcessed() throws Exception {
        AnalyzeBAMSandAllenGenes.Direction direction = AnalyzeBAMSandAllenGenes.Direction.ANYDIRECTION;
        boolean removeNonExp = true;
        boolean useVirtual = true;
        ABAMSDataMatrix geneRegionExpressionMatrix = ExpressionMatrixPairFactory.connectivityAndExpression( direction,
                useVirtual, removeNonExp ).getMatrixB();
        Util.writeRTable( SetupParameters.getDataFolder() + "ABAMSProcessedMatrix.txt", geneRegionExpressionMatrix );
    }

    /**
     * @param args
     */
    public static void main( String[] args ) throws Exception {
        // getAll();
        // getAllenProcessed();
        getABAMSProcessed();
        System.exit( 1 );
        AnalyzeBAMSandAllenGenes.Direction direction = AnalyzeBAMSandAllenGenes.Direction.ANYDIRECTION;
        boolean removeNonExp = true;
        ABAMSDataMatrix geneRegionExpressionMatrix = ExpressionMatrixPairFactory.getEnergyMatrix( direction,
                removeNonExp );

        // loads in the brain region heirarchy
        StructureCatalogLoader loader = new StructureCatalogLoader();
        Collection<String> nonLeafs = ( Collection<String> ) Util.subtract( geneRegionExpressionMatrix.getColNames(),
                loader.getLeafs() );
        geneRegionExpressionMatrix = geneRegionExpressionMatrix.removeColumns( nonLeafs );

        GeneFilter filter;
        List<String> removeRows;

        // Coronal only?
        // GeneFilter filter = new PlaneRemoveFilter( PlaneRemoveFilter.Plane.SAGITTAL );
        // List<String> removeRows = filter.getRowsToRemove( geneRegionExpressionMatrix );
        // log.info( "Removed:" + removeRows.size() + " Name:" + filter.getName() );
        // geneRegionExpressionMatrix = geneRegionExpressionMatrix.removeRows( removeRows );

        filter = new NaNGeneFilter();
        removeRows = filter.getRowsToRemove( geneRegionExpressionMatrix );
        log.info( "Removed:" + removeRows.size() + " Name:" + filter.getName() );
        geneRegionExpressionMatrix = geneRegionExpressionMatrix.removeRows( removeRows );

        System.out.println( "Rows(Gene assays):" + geneRegionExpressionMatrix.rows() + " Cols(Brain regions):"
                + geneRegionExpressionMatrix.columns() );

        Util.writeRTable( "/grp/java/workspace/BAMSandAllen/data/AllenDataForWyethLab."
                + geneRegionExpressionMatrix.columns() + "regions.x." + geneRegionExpressionMatrix.rows()
                + ".assays.txt", geneRegionExpressionMatrix );
        System.exit( 1 );
        // write out
        log.info( geneRegionExpressionMatrix.getColNames() );
        // iterate genes
        for ( String geneAssay : geneRegionExpressionMatrix.getRowNames() ) {
            if ( geneAssay.startsWith( "Rbpms2" ) ) {
                log.info( geneAssay );
                log.info( geneRegionExpressionMatrix.getByKeys( geneAssay, "Facial motor nucleus" ) );
                double[] expression = geneRegionExpressionMatrix.getRowByName( geneAssay );
                // iterage regions
                for ( String brainregion : geneRegionExpressionMatrix.getColNames() ) {
                    geneRegionExpressionMatrix.getByKeys( geneAssay, brainregion );
                }
            }
        }
    }
}
