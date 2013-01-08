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

import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.ABAMSDataMatrix;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes;
import ubic.BAMSandAllen.Util;
import ubic.BAMSandAllen.AllenDataLoaders.StructureCatalogLoader;
import ubic.BAMSandAllen.geneFilters.GeneFilter;
import ubic.BAMSandAllen.geneFilters.PlaneRemoveFilter;
import ubic.BAMSandAllen.optimize.GreedyMultiThreaded;

public class ExpressionLiteratureMatrixPair extends LiteratureMatrixPair {
    private static Log log = LogFactory.getLog( ExpressionLiteratureMatrixPair.class.getName() );

    public ExpressionLiteratureMatrixPair() throws Exception {
        super( Nomenclature.ABA );

        AnalyzeBAMSandAllenGenes.Direction direction = AnalyzeBAMSandAllenGenes.Direction.INCOMING;
        boolean removeNonExp = true;
        ABAMSDataMatrix geneRegionExpressionMatrix = ExpressionMatrixPairFactory.getEnergyMatrix( direction,
                removeNonExp );

        // loads in the brain region heirarchy
        StructureCatalogLoader loader = new StructureCatalogLoader();
        Collection<String> nonLeafs = ( Collection<String> ) Util.subtract( geneRegionExpressionMatrix.getColNames(),
                loader.getLeafs() );
        matrixA = geneRegionExpressionMatrix.removeColumns( nonLeafs );

        log.info( "Rows(Gene assays):" + geneRegionExpressionMatrix.rows() + " Cols(Brain regions):"
                + geneRegionExpressionMatrix.columns() );
    }

    public void removeSaggitalGenes() {
        GeneFilter filter = new PlaneRemoveFilter( PlaneRemoveFilter.Plane.SAGITTAL );
        List<String> removeRows = filter.getRowsToRemove( matrixA );
        log.info( "Removed:" + removeRows.size() + " Name:" + filter.getName() );
        matrixA = matrixA.removeRows( removeRows );
    }

    public static void main( String[] args ) throws Exception {
        ExpressionLiteratureMatrixPair pair = new ExpressionLiteratureMatrixPair();
        pair.printDimensions();
        pair.run();
        log.info( pair.getMatrixB().getRowByName( "the" ) );
        Set<String> rows = new HashSet<String>( pair.getMatrixB().getRowNames() );
        pair.removeLiteratureOnes();
        pair.printDimensions();
        log.info( "subtracted:" + Util.subtract( rows, pair.getMatrixB().getRowNames() ) );
        // pair.writeImages();
        pair.removeSaggitalGenes();
        //pair.switchMatrices();

        //System.exit( 1 );
        int threads = 16;
        GreedyMultiThreaded remover = new GreedyMultiThreaded( pair, threads );
        boolean keepSign = true; // should it increase or decrease the correlation?
        remover.run( pair.getMatrixB().rows(), keepSign );
    }

}
