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

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes.Direction;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenExpressionMatrixPair;
import ubic.BAMSandAllen.MatrixPairs.ExpressionMatrixPairFactory;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenPartialExpressionMatrixPair.RegressMatrix;
import ubic.BAMSandAllen.geneFilters.AGEAGeneFilter;

public class ForAGEAGenes {
    private static Log log = LogFactory.getLog( ForAGEAGenes.class.getName() );

    /**
     * @param args
     */
    public static void main( String[] args ) throws Exception {
        Direction direction;
        direction = AnalyzeBAMSandAllenGenes.Direction.OUTGOING;
        boolean slow = false;
        boolean logDistance = true;
        RegressMatrix regressType = RegressMatrix.BOTH;
        boolean useVirtual = true;
        boolean removeNonExp = true;

        ConnectivityAndAllenExpressionMatrixPair partialMantelPair = ExpressionMatrixPairFactory.connectivityPartial(
                direction, slow, regressType, useVirtual, removeNonExp, logDistance );
        boolean imageSeries = true;
        partialMantelPair.applyGeneFilter( new AGEAGeneFilter( imageSeries ) );
        log.info( partialMantelPair.getCorrelation() );

    }

}
