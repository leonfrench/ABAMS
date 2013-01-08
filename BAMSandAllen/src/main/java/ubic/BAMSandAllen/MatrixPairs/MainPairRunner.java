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

import java.util.Map;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes;
import ubic.BAMSandAllen.SetupParameters;
import ubic.BAMSandAllen.Util;
import ubic.BAMSandAllen.AllenDataLoaders.AllenCSVLoader;
import ubic.BAMSandAllen.ClassSelectors.BrainRegionClassSelector;
import ubic.basecode.dataStructure.params.ParamKeeper;
import ubic.basecode.dataStructure.params.ParameterGrabber;

public class MainPairRunner {
    private static Log log = LogFactory.getLog( MainPairRunner.class.getName() );

    public static void main( String[] args ) throws Exception {

        // int iterations = 1000;

        int iterations = 1000;
        double[] zeroReplacements = { Double.NaN };
        // new energy low= 0.0000000000647
        boolean[] doLogs = { true };
        boolean[] do1pLogs = { false };
        boolean[] squares = { false };
        // String[] matrices = { "NewEnergies", "NewLevels", "NewDensity", "OriginalEnergies", "OriginalLevels",
        // "OriginalDensity", "MinerLevel", "MinerMean" };
        // String[] matrices = { "NewEnergies", "NewLevels", "NewDensity", "OriginalEnergies", "MinerLevel", "MinerMean"
        // };
        String[] matrices = { "NewEnergies" };

        // // config that gives bad p-value
        // double[] zeroReplacements = { Double.NaN, 0, 0.0000000000647 };
        // boolean[] doLogs = { true };
        // boolean[] do1pLogs = { false };
        // boolean[] squares = { true };
        // String[] matrices = { "NewEnergies" };

        // //NewEnergies NewLevels NewDensity

        // double[] zeroReplacements = { 0, 0.01d, 0.001d, Double.NaN };
        // boolean[] doLogs = { true, false };
        // boolean[] squares = { true, false };
        // String[] matrices = { "MinerLevel", "MinerMean", "OriginalEnergies", "NewEnergies" };

        ParamKeeper result = new ParamKeeper();

        for ( double zeroReplacement : zeroReplacements ) {
            for ( boolean doLog : doLogs ) {
                for ( boolean do1pLog : do1pLogs ) {
                    for ( boolean square : squares ) {
                        for ( String matrixname : matrices ) {
                            ConnectivityAndAllenExpressionMatrixPair matrixPair = new ConnectivityAndAllenExpressionMatrixPair(
                                    new BrainRegionClassSelector(), doLog, do1pLog, square, zeroReplacement, matrixname, AnalyzeBAMSandAllenGenes.Direction.OUTGOING );
                            double correlation = matrixPair.run();
                            double pvalue = matrixPair.test( iterations );
                            double zeroes = Util.countValues( matrixPair.matrixB, 0d );
                            double NaNs = Util.countValues( matrixPair.matrixB, Double.NaN );

                            Map<String, String> params = ParameterGrabber.getParams(
                                    ConnectivityAndAllenExpressionMatrixPair.class, matrixPair );
                            params.put( "iterations", iterations + "" );
                            params.put( "pvalue", pvalue + "" );
                            params.put( "correlation", correlation + "" );
                            params.put( "zeroes", zeroes + "" );
                            params.put( "NaNs", NaNs + "" );
                            params.put( "rows", matrixPair.matrixB.rows() + "" );
                            params.put( "cols", matrixPair.matrixB.columns() + "" );
                            params.put( "connections", matrixPair.getConnectionCount() + "" );
                            Set<String> genes = AllenCSVLoader.getAllGenes( matrixPair.matrixB.getRowNames() );
                            params.put( "genes", genes.size() + "" );
                            // matrixPair.writeRMatrices();
                            result.addParamInstance( params );
                        }
                    }
                }
                System.out.println( "So far:" );
                System.out.println( result.toCSVString() );
                // writeExcel
            }
        }
        // System.out.println( result.toCSVString() );
        result.writeExcel( SetupParameters.getDataFolder() + "MainPairRunner." + System.currentTimeMillis() + ".xls" );
        log.info( "Wrote to excel file" );
    }
}
