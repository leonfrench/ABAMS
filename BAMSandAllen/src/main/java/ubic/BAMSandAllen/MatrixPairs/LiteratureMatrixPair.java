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
import java.util.HashSet;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.ABAMSDataMatrix;
import ubic.BAMSandAllen.SetupParameters;
import ubic.BAMSandAllen.adjacency.CorrelationAdjacency;
import ubic.BAMSandAllen.geneFilters.RemoveAllOneRows;
import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;

public abstract class LiteratureMatrixPair extends MatrixPair {
    private static Log log = LogFactory.getLog( ConnectivityLiteratureMatrixPair.class.getName() );

    public enum Nomenclature {
        ABA, BAMS
    };

    public String getRegionsForWord( String word ) {
        int count = 0;
        String result = "";
        String notIn = "";
        for ( String region : matrixB.getColNames() ) {
            double value = matrixB.getByKeys( word, region );
            if ( value == 1d ) {
                count++;
                result += region + ", ";
            } else
                notIn += region + ", ";
        }
        result = word + ": " + count + " regions (" + result + ")";
        if ( ( matrixB.columns() - count ) < 15 ) result += " not in:" + notIn;
        return result;
    }

    Nomenclature nomenclature;

    public LiteratureMatrixPair( Nomenclature nomenclature ) throws Exception {
        this.nomenclature = nomenclature;
        isInSameSpace = false;

        String filename = SetupParameters.config.getString( "abams.dataFolder" );
        filename += nomenclature.toString() + " concept.matrix.cache";

        String nameB = "Literature" + nomenclature.toString();

        ObjectInputStream in = new ObjectInputStream( new FileInputStream( filename ) );
        DoubleMatrix<String, String> matrixBRead = ( DenseDoubleMatrix<String, String> ) in.readObject();
        in.close();

        // adjacency methods, needs to be paramertized/chosen
        matrixB = new ABAMSDataMatrix( matrixBRead, nameB, new CorrelationAdjacency( matrixBRead ) );
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

    public void removeLiteratureOnes() {
        // RemoveLowVarianceGeneFilter filter = new RemoveLowVarianceGeneFilter( 0 );
        RemoveAllOneRows filter = new RemoveAllOneRows();
        matrixB = matrixB.removeRows( filter.getRowsToRemove( matrixB ) );
    }

    public double run() throws Exception {
        slimMatrices();
        sameSpace();
        double correl = getCorrelation( false );
        printDimensions();
        log.info( matrixA.getName() + " and " + matrixB.getName() + " correlation is " + correl );
        test( 1000, false );
        // writeImages();
        return correl;
    }

}
