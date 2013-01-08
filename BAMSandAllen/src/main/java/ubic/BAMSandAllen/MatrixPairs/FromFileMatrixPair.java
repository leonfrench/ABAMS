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

import java.io.File;
import java.util.HashSet;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.ABAMSDataMatrix;
import ubic.BAMSandAllen.SetupParameters;
import ubic.BAMSandAllen.adjacency.AdjacencyCompute;
import ubic.BAMSandAllen.adjacency.CorrelationAdjacency;
import ubic.BAMSandAllen.adjacency.DotProductAdjacency;
import ubic.BAMSandAllen.adjacency.ManhattanAdjacency;
import ubic.BAMSandAllen.optimize.GreedyMultiThreaded;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.io.reader.DoubleMatrixReader;

public class FromFileMatrixPair extends MatrixPair {
    private static Log log = LogFactory.getLog( FromFileMatrixPair.class.getName() );

    public FromFileMatrixPair( String fileNameA, String fileNameB ) throws Exception {
        super();
        isInSameSpace = true;
        DoubleMatrixReader reader = new DoubleMatrixReader();
//        reader.setTopLeft( false );

        DoubleMatrix<String, String> matrixARead = reader.read( fileNameA );
        DoubleMatrix<String, String> matrixBRead = reader.read( fileNameB );

        String nameA = new File( fileNameA ).getName();
        String nameB = new File( fileNameB ).getName();

        // adjacency methods, needs to be paramertized
        this.matrixA = new ABAMSDataMatrix( matrixARead, nameA, new DotProductAdjacency() );
        this.matrixB = new ABAMSDataMatrix( matrixBRead, nameB, new CorrelationAdjacency( matrixBRead ) );
    }

    public Set<String> convertANametoB( String aName ) {
        Set<String> result = new HashSet<String>();
        result.add( aName );
        return result;
    }

    public Set<String> convertBNametoA( String bName ) {
        Set<String> result = new HashSet<String>();
        result.add( bName );
        return result;
    }

    public double run() throws Exception {
        sameSpace();
        double correl = getCorrelation( false );
        printDimensions();
        log.info( matrixA.getName() + " and " + matrixB.getName() + " correlation is " + correl );
        testBoth( 1000, false );
        writeImages();
        return correl;
    }

    /**
     * @param args
     */
    public static void main( String[] args ) throws Exception {

        String fileA = SetupParameters.getDataFolder() + "mouse.ann.mat";
        String fileB = SetupParameters.getDataFolder() + "mouse.rank.mat";

        MatrixPair pair = new FromFileMatrixPair( fileA, fileB );
        pair.run();

        int threads = 8;
        GreedyMultiThreaded remover = new GreedyMultiThreaded( pair, threads );
        boolean keepSign = true; // should it increase or decrease the correlation?
        remover.run( 6882, keepSign );
    }

}
