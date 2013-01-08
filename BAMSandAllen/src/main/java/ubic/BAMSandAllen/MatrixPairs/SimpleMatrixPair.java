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
import ubic.BAMSandAllen.BAMSDataLoaders.BAMSFromNIFDataLoader;

public class SimpleMatrixPair extends MatrixPair {
    protected static Log log = LogFactory.getLog( SimpleMatrixPair.class );

    public SimpleMatrixPair( ABAMSDataMatrix a, ABAMSDataMatrix b ) {
        matrixA = a;
        matrixB = b;
        isInSameSpace = false;
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
        slimMatrices();
        sameSpace();
        double correl = getCorrelation( false );
        printDimensions();
        log.info( matrixA.getName() + " and " + matrixB.getName() + " correlation is " + correl );
        test( 1000, false );
        // writeImages();
        return correl;
    }

    /**
     * @param args
     */
    public static void main( String[] args ) {
        // TODO Auto-generated method stub

    }

}
