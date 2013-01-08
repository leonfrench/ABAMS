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
import ubic.BAMSandAllen.Util;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes.Direction;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenExpressionMatrixPair;
import ubic.BAMSandAllen.MatrixPairs.ExpressionMatrixPairFactory;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenPartialExpressionMatrixPair.RegressMatrix;
import ubic.BAMSandAllen.geneFilters.EphrinGeneFilter;
import ubic.BAMSandAllen.geneFilters.GeneFilter;

public class EstrogenTests {
    private static Log log = LogFactory.getLog( EstrogenTests.class.getName() );

    GeneFilter f;

    public EstrogenTests( GeneFilter f ) {
        this.f = f;
    }

    public static void main( String[] args ) throws Exception {
        String toPrint = "";
        EstrogenTests e;

        // e = new EstrogenTests( new PrefixGeneFilter("Sema") );
        // toPrint += e.run();

        // e = new EstrogenTests( null );
        // toPrint += e.run();

        e = new EstrogenTests( new EphrinGeneFilter() );
        toPrint += e.run();

        // e = new EstrogenTests( new EstrogenGeneFilter() );
        // toPrint += e.run();
        System.out.println( toPrint );
    }

    /**
     * @param args
     */
    public String run() throws Exception {
        String result = "\n";
        if ( f != null )
            result += f.getName() + "\n";
        else
            result += "No gene Filter\n";

        boolean removeNonExp = true;
        boolean useVirtual = true;

        Direction direction = AnalyzeBAMSandAllenGenes.Direction.OUTGOING;
        ConnectivityAndAllenExpressionMatrixPair pair = ExpressionMatrixPairFactory.connectivityAndExpression(
                direction, useVirtual, removeNonExp );
        result += "Outgoing:\n";
        result += test( pair );

        // ///////////////////
        // pair.writeImages();
        // pair.writeRMatrices();
        // System.exit( 1 );
        // ///////////////////

        result += "Incoming:\n";
        direction = AnalyzeBAMSandAllenGenes.Direction.INCOMING;
        pair = ExpressionMatrixPairFactory.connectivityAndExpression( direction, useVirtual, removeNonExp );
        result += test( pair );

        log.info( "Number of genes:" + Util.getUniqueGenes( pair.getMatrixBDataRows() ) );

        boolean slow = false;
        boolean logDistance = true;
        RegressMatrix regressType = RegressMatrix.BOTH;

        result += "Outgoing Partial:\n";
        direction = AnalyzeBAMSandAllenGenes.Direction.OUTGOING;
        pair = ExpressionMatrixPairFactory.connectivityPartial( direction, slow, regressType, useVirtual, removeNonExp,
                logDistance );
        result += test( pair );

        result += "Incoming Partial:\n";
        direction = AnalyzeBAMSandAllenGenes.Direction.INCOMING;
        pair = ExpressionMatrixPairFactory.connectivityPartial( direction, slow, regressType, useVirtual, removeNonExp,
                logDistance );
        result += test( pair );

        boolean run = false;
        boolean square = true;
        useVirtual = false;
        direction = AnalyzeBAMSandAllenGenes.Direction.OUTGOING;
        pair = ExpressionMatrixPairFactory.connectivityAndExpression( direction, useVirtual, removeNonExp, run, square );

        result += "Outgoing Direct:\n";
        pair.makeDirect();
        pair.run();
        result += test( pair );

        result += "Incoming Direct:\n";
        direction = AnalyzeBAMSandAllenGenes.Direction.INCOMING;
        pair = ExpressionMatrixPairFactory.connectivityAndExpression( direction, useVirtual, removeNonExp, run, square );
        pair.makeDirect();
        pair.run();
        result += test( pair );

        result += "Number of image series sets:" + pair.getMatrixBDataRows().size() + "\n";
        result += "Number of genes:" + Util.getUniqueGenes( pair.getMatrixBDataRows() ).size() + "\n";
        if ( f != null ) result += "Genes:" + Util.getUniqueGenes( pair.getMatrixBDataRows() ) + "\n";

        // do direct
        return result;
        // direct controlled for space - ugh

    }

    public String test( ConnectivityAndAllenExpressionMatrixPair pair ) throws Exception {
        String result = "";
        if ( f != null ) pair.applyGeneFilter( f );
        result += "Correlation:" + pair.getCorrelation( true ) + " ";
        result += "P-value:" + pair.test( 1000, false ) + "\n";
        return result;
    }

}
