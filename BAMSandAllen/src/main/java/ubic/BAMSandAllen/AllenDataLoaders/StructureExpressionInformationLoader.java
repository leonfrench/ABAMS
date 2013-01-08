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

import java.io.FileReader;
import java.io.FileWriter;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.io.writer.MatrixWriter;
import au.com.bytecode.opencsv.CSVReader;

public class StructureExpressionInformationLoader extends AllenCSVLoader {
    private static Log log = LogFactory.getLog( StructureExpressionInformationLoader.class.getName() );

    // cvs header information from BrainFeatureTable.csv
    public static final int geneNamePosition = 0;
    public static final int imageIDPosition = 1;
    public static final int densityPosition = 5;
    public static final int levelPosition = 6;
    public static final int energyPosition = 7;

    public StructureExpressionInformationLoader() throws Exception {
        super();
        headerLine = true;
        sep = ',';
        filename = "/home/leon/temp/allen/BrainFeatureTable.csv";
        regionNamePosition = 4;
        init();
    }

    public String getRowName( String line[] ) {
        return line[geneNamePosition] + "[" + line[imageIDPosition] + "]";
    }

    public static void main( String[] args ) throws Exception {
        // AllenCSVLoader loader = new NewStructureExpressionInformationLoader();
        AllenCSVLoader loader = new AllenMinerLoader();
        // loader.convertToRTable( "tempEnergyMatrix.txt", energyPosition );
        // loader.convertToRTable( "tempLevelMatrix.txt", levelPosition );
        // loader.convertToRTable( "tempDensityMatrix.txt", densityPosition );
        List<String> ubiq = loader.getGeneListFromFile( "/grp/java/workspace/BAMSandAllen/data/ABAUbiquitous.txt" );
        List<String> nonExp = loader.getGeneListFromFile( "/grp/java/workspace/BAMSandAllen/data/ABANonexpressed.txt" );
        System.out.println( ubiq.size() );
        System.out.println( nonExp.size() );
        Set<String> genes2Keep = loader.imageSeriesNames;
        log.info( "Before ubiq and nonExp:" + loader.imageSeriesNames.size() );
        genes2Keep.removeAll( ubiq );
        genes2Keep.removeAll( nonExp );
        log.info( "After:" + loader.imageSeriesNames.size() );
        log.info( "After genes:" + loader.getAllGenes().size() );
        log.info( genes2Keep.size() );

        // loader.convertToRTable( "tempEnergyMatrixReduced.txt", energyPosition, new LinkedList<String>( genes2Keep )
        // );
    }
}
