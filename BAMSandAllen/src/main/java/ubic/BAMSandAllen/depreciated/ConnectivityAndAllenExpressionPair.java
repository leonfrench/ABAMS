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
package ubic.BAMSandAllen.depreciated;

import java.io.FileReader;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import org.apache.commons.lang.time.StopWatch;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import au.com.bytecode.opencsv.CSVReader;

import ubic.BAMSandAllen.SetupParameters;
import ubic.BAMSandAllen.Util;
import ubic.BAMSandAllen.AllenDataLoaders.AllenCatalogMatrices;
import ubic.BAMSandAllen.ClassSelectors.BrainRegionClassSelector;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenDataPair;

@Deprecated
public class ConnectivityAndAllenExpressionPair { //extends ConnectivityAndAllenDataPair {
//    private static Log log = LogFactory.getLog( ConnectivityAndAllenExpressionPair.class.getName() );
//
//    public ConnectivityAndAllenExpressionPair( BrainRegionClassSelector selector, boolean doLog, boolean square,
//            double zeroReplacement ) throws Exception {
//        super( selector, square );
//
//        AllenCatalogMatrices allenMatrices = new AllenCatalogMatrices();
//        matrixB = allenMatrices.getEnergies();
//
//        processExpression();
//
//        log.info( "got Matrix B - log Expression" );
//
//        if ( doLog ) {
//            matrixB = Util.logMatrix( matrixB, zeroReplacement );
//        }
//    }
//
//    protected void processExpression() {
//        // any rows with zero exp?
//
//        // remove cols with no expression
//        matrixB = removeColumns( matrixB, Util.findZeroColumns( matrixB ) );
//
//        try {
//            List<String> ubiq = getGeneRowNamesFromFile( SetupParameters.getDataFolder() + "ABAUbiquitous.txt" );
//            List<String> nonExp = getGeneRowNamesFromFile( SetupParameters.getDataFolder() + "ABANonexpressed.txt" );
//            log.info( "Ubiquitous:" + ubiq.size() );
//            log.info( "Non expressors:" + nonExp.size() );
//            matrixB = Util.removeRows( matrixB, new HashSet<String>( nonExp ) );
//            matrixB = Util.removeRows( matrixB, new HashSet<String>( ubiq ) );
//        } catch ( Exception e ) {
//            e.printStackTrace();
//            log.info( "Failed to load ubiq and nonexp genes" );
//        }
//    }
//
//    public List<String> getGeneRowNamesFromFile( String filename ) throws Exception {
//        List<String> genes = new LinkedList<String>();
//        CSVReader reader = new CSVReader( new FileReader( filename ) );
//        String[] line = reader.readNext();
//        while ( ( line = reader.readNext() ) != null ) {
//            genes.add( line[0] );
//        }
//        reader.close();
//        List<String> resultGenes = new LinkedList<String>();
//        for ( String gene : genes ) {
//            for ( String row : matrixB.getRowNames() ) {
//                if ( row.startsWith( gene + "[" ) ) {
//                    resultGenes.add( row );
//                }
//            }
//        }
//        return resultGenes;
//    }
//
//    public static void main( String[] args ) throws Exception {
//        double zeroReplacement = 0;
//        boolean useSquareConnectivity = false;
//
//        ConnectivityAndAllenExpressionPair x = new ConnectivityAndAllenExpressionPair( new BrainRegionClassSelector(),
//                true, useSquareConnectivity, zeroReplacement );
//        x.run();
//        x.writeRMatrices();
//        // System.exit( 1 );
//        // x.writeImages();
//        StopWatch watch = new StopWatch();
//        watch.start();
//        x.test( 100 );
//        log.info( watch.toString() );
//
//        useSquareConnectivity = !useSquareConnectivity;
//        x = new ConnectivityAndAllenExpressionPair( new BrainRegionClassSelector(), true, useSquareConnectivity,
//                zeroReplacement );
//        x.run();
//        // x.writeRMatrices();
//        // x.writeImages();
//        watch = new StopWatch();
//        watch.start();
//        x.test( 1000 );
//        log.info( watch.toString() );
//
//    }
}
