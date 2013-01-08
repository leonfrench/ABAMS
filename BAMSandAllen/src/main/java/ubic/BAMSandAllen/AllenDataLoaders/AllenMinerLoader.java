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
package ubic.BAMSandAllen.AllenDataLoaders;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.ABAMSDataMatrix;
import ubic.BAMSandAllen.adjacency.CorrelationAdjacency;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.math.MatrixStats;

public class AllenMinerLoader extends AllenCSVLoader {
    private static Log log = LogFactory.getLog( StructureExpressionInformationLoader.class.getName() );

    public static final int geneNamePosition = 0;
    public static final int slicePosition = 1;
    public static final int xprFilePosition = 2;
    public static final int numROIPointsPosition = 4;
    public static final int numExpressingROIPointsPosition = 5;
    public static final int totalExpressionInROI = 6;
    public static final int meanExpressionPosition = 7;

    public AllenMinerLoader() throws Exception {
        super();
        headerLine = false;
        filename = "/home/leon/temp/allen/allenminer/manuscript_data/ABA_region_expression/all_regions.roi_list.out";
        filename = "/home/leon/temp/allen/allenminer/manuscript_data/ABA_region_expression/all_regions.roi_list.20100123.out";

        regionNamePosition = 3;
        sep = '\t';

        init();
    }

    public static void main( String[] args ) throws Exception {
        AllenMinerLoader miner = new AllenMinerLoader();
        log.info( miner.getRegionNames() );
        DoubleMatrix<String, String> matrix = miner.getMatrix( AllenMinerLoader.totalExpressionInROI );
        ABAMSDataMatrix aMatrix = new ABAMSDataMatrix( matrix, "AllenMiner", new CorrelationAdjacency( matrix ) );
        boolean[][] nans = MatrixStats.nanStatusMatrix( aMatrix.getRawMatrix() );
        int nanCount = 0;
        int total = 0;
        for ( int i = 0; i < nans.length; i++ ) {
            for ( int j = 0; j < nans[0].length; j++ ) {
                if ( nans[i][j] ) nanCount++;
                total++;
            }
        }
        log.info( "Total:" + total );
        log.info( "NaNs:" + nanCount );
        // allenMatrices.getMinerEnergies();
    }

    public String getRowName( String[] line ) {
        return line[geneNamePosition] + "[" + line[xprFilePosition] + "]." + line[slicePosition];
    }
}
