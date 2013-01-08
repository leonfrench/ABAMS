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

import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes;
import ubic.BAMSandAllen.SetupParameters;
import ubic.BAMSandAllen.Util;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes.Direction;
import ubic.BAMSandAllen.MatrixPairs.AllenMatrixPair;
import ubic.BAMSandAllen.MatrixPairs.AllenMatrixPairFactory;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenExpressionMatrixPair;
import ubic.BAMSandAllen.MatrixPairs.ExpressionMatrixPairFactory;
import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;

public class CreateAllenDegreeVector {
    private static Log log = LogFactory.getLog( CreateAllenDegreeVector.class.getName() );

    /**
     * @param args
     */
    public static void main( String[] args ) throws Exception {
        Direction direction = AnalyzeBAMSandAllenGenes.Direction.APPENDED;

        boolean logDistance = true;
        boolean useVirtual = true;
        boolean run = true;
        boolean square = false;
        boolean removeNonExp = true;

        AllenMatrixPair spacePair = AllenMatrixPairFactory.getSpaceAndExpressionPair( direction, removeNonExp,
                logDistance );
        spacePair.printDimensions();
        spacePair.writeRMatrices();
        System.exit(1);

        ConnectivityAndAllenExpressionMatrixPair pair = ExpressionMatrixPairFactory.connectivityAndExpression(
                direction, useVirtual, removeNonExp, run, square );

        boolean removeNan = true;
        DoubleMatrix<String, String> degrees = Util.columnSums( pair.getMatrixA(), removeNan );

        List<String> regionNamesBAMS = degrees.getColNames();

        DoubleMatrix<String, String> result = new DenseDoubleMatrix<String, String>( 1, spacePair.getMatrixA()
                .getColNames().size() );
        result.setColumnNames( new LinkedList<String>( spacePair.getMatrixA().getColNames() ) );
        result.addRowName( "degree" );

        Set<String> regionNamesAllen = new HashSet<String>();
        for ( String BAMSName : regionNamesBAMS ) {
            Set<String> allenMapped = pair.convertANametoB( BAMSName );
            if ( allenMapped.size() > 1 ) {
                throw new RuntimeException( BAMSName + " Has more than one Allen mapping" );
            }
            if ( allenMapped.size() == 0 ) {
                throw new RuntimeException( BAMSName + " Has zero Allen mapping!" );
            }
            String allenRegionName = allenMapped.iterator().next();
            regionNamesAllen.addAll( allenMapped );
            result.setByKeys( "degree", allenRegionName, degrees.getByKeys( "Sums", BAMSName ) );
        }
        log.info( "Total Allen regions:" + regionNamesAllen.size() );

        Util.writeRTable( SetupParameters.getDataFolder() + "AllenDegrees" + direction.toString() + ".txt", result );

    }
}
