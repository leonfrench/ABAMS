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
package ubic.BAMSandAllen.FocusedAnalysis;

import ubic.BAMSandAllen.ABAMSDataMatrix;
import ubic.BAMSandAllen.Util;
import ubic.BAMSandAllen.AllenDataLoaders.AllenCatalogMatrices2;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes.Direction;
import ubic.BAMSandAllen.MatrixPairs.ExpressionMatrixPairFactory;
import ubic.BAMSandAllen.adjacency.CorrelationAdjacency;
import ubic.BAMSandAllen.geneFilters.NonExpFilter;
import ubic.BAMSandAllen.geneFilters.PlaneRemoveFilter;
import ubic.BAMSandAllen.geneFilters.PlaneRemoveFilter.Plane;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.math.CorrelationStats;

public class SimpleBrainRegionMaskTest {

    public static void main( String args[] ) throws Exception {
        // load expression data
        String matrixName = "NewEnergies";

        AllenCatalogMatrices2 allenMatrices = new AllenCatalogMatrices2();
        DoubleMatrix<String, String> dataMatrix = allenMatrices.getFromDisk( matrixName );
        ABAMSDataMatrix energy = new ABAMSDataMatrix( dataMatrix, matrixName, new CorrelationAdjacency( dataMatrix ) );
        for ( String name : energy.getColNames() ) {
//            System.out.println( name );
        }

        double[] mask = energy.getRow( 1 );
        for ( int i = 0; i < mask.length; i++ ) {
            mask[i] = 0;
        }
        int col;
        col = energy.getColIndexByName( "Spinal nucleus of the trigeminal_ caudal part" );
        mask[col] = 1;
        col = energy.getColIndexByName( "Spinal nucleus of the trigeminal_ interpolar part" );
        mask[col] = 1;
        col = energy.getColIndexByName( "Spinal nucleus of the trigeminal_ oral part" );
        mask[col] = 1;
        col = energy.getColIndexByName( "Ventral posterior complex of the thalamus" );
        mask[col] = 1;
        col = energy.getColIndexByName( "Periaqueductal gray" );
        mask[col] = 1;
        col = energy.getColIndexByName( "Principal sensory nucleus of the trigeminal" );
        mask[col] = 1;
        
//        energy = energy.applyRowFilter( new PlaneRemoveFilter(Plane.SAGITTAL) );
        energy = energy.applyRowFilter( new NonExpFilter() );

        for ( String gene : energy.getRowNames() ) {
            double[] exp = energy.getRowByName( gene );
            double pearson = CorrelationStats.correl( mask, exp );
            double spearman = Util.spearmanCorrel( mask, exp );
            System.out.println( gene + "," + spearman + "," + pearson + ","
                    + CorrelationStats.pvalue( pearson, exp.length ) );

        }
        // Spinal nucleus of the trigeminal_ caudal part
        // Spinal nucleus of the trigeminal_ interpolar part
        // Spinal nucleus of the trigeminal_ oral part
        // Ventral posterior complex of the thalamus
        // Periaqueductal gray
        // Principal sensory nucleus of the trigeminal

    }
    // public double getDegreeCorrelation( String rowName ) {
    // double[] exp = matrixB.getRowByName( rowName );
    // DoubleMatrix<String, String> connectionDegrees = Util.columnSums( matrixA );
    // double[] degrees = connectionDegrees.getRow( 0 );
    // return CorrelationStats.correl( degrees, exp );
    // }
}
