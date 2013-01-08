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
import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.StringTokenizer;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.Util;
import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import au.com.bytecode.opencsv.CSVReader;

public abstract class AllenCSVLoader {
    private static Log log = LogFactory.getLog( AllenCSVLoader.class.getName() );

    boolean headerLine;
    char sep;
    public int regionNamePosition = -1;
    String filename;
    StructureCatalogLoader nameResolver;
    Set<String> imageSeriesNames;
    Set<String> regionNames;

    public Set<String> getRegionNames() {
        return regionNames;
    }

    public AllenCSVLoader() throws Exception {
        imageSeriesNames = new HashSet<String>();
        regionNames = new HashSet<String>();
        nameResolver = new StructureCatalogLoader();
    }

    protected String getFullRegionName( String[] line ) {
        String acro = line[regionNamePosition];
        if ( acro.equals( "internal_background" ) ) {
            acro = "Brain";
        }
        // return "\"" + nameResolver.getFullName( acro ) + "\"";
        return nameResolver.getFullName( acro );
    }

    public void convertToRTable( String filename, int columnPosition ) throws Exception {
        DoubleMatrix<String, String> matrix = getMatrix( columnPosition );
        Util.writeRTable( filename, matrix );
    }

    public void convertToRTable( String filename, int columnPosition, List<String> geneList ) throws Exception {
        DoubleMatrix<String, String> matrix = getMatrix( columnPosition );
        Util.writeRTable( filename, matrix );
    }

    public List<String> getGeneListFromFile( String filename ) throws Exception {
        List<String> genes = new LinkedList<String>();
        CSVReader reader = new CSVReader( new FileReader( filename ) );
        String[] line = reader.readNext();
        while ( ( line = reader.readNext() ) != null ) {
            // System.out.println( line[0] );
            genes.add( line[0] );
        }
        reader.close();
        List<String> resultGenes = new LinkedList<String>();
        for ( String gene : genes ) {
            for ( String row : imageSeriesNames ) {
                if ( row.startsWith( gene + "[" ) ) {
                    resultGenes.add( row );
                }
            }
        }
        return resultGenes;
    }

    public DoubleMatrix<String, String> getMatrix( int columnPosition ) throws Exception {
        // we want all the genes
        return getMatrix( columnPosition, new LinkedList<String>( imageSeriesNames ) );
    }

    public DoubleMatrix<String, String> getMatrix( int columnPosition, List<String> genes2Keep ) throws Exception {
        DoubleMatrix<String, String> resultMatrix = new DenseDoubleMatrix<String, String>( imageSeriesNames.size(),
                regionNames.size() );
        resultMatrix.setRowNames( genes2Keep );
        resultMatrix.setColumnNames( new LinkedList<String>( regionNames ) );

        CSVReader reader = new CSVReader( new FileReader( filename ), sep );
        // eat the header line
        String[] line;
        if ( headerLine ) reader.readNext();
        while ( ( line = reader.readNext() ) != null ) {
            String rowName = getRowName( line );
            if ( genes2Keep.contains( rowName ) ) {
                // convert the acronym to the full name
                String colName = getFullRegionName( line );
                // quote it
                // the actuall value we want in the matrix
                String output = line[columnPosition];
                Double value;
                if ( output.equals( "NA" ) ) {
                     value = Double.NaN;
                    //value = 0d;
                } else
                    value = Double.parseDouble( output );

                resultMatrix.setByKeys( rowName, colName, value );
            }
        }
        reader.close();
        return resultMatrix;
    }

    public void init() throws Exception {
        CSVReader reader = new CSVReader( new FileReader( filename ), sep );
        String[] line;
        // eat the header line
        int linecount = 0;
        if ( headerLine ) {
            linecount++;
            line = reader.readNext();
        }
        while ( ( line = reader.readNext() ) != null ) {
            linecount++;
            imageSeriesNames.add( getRowName( line ) );
            // log.info( line.length );
            // log.info( line[regionNamePosition] );
            // log.info( getQuotedRegionName( line ) );
            regionNames.add( getFullRegionName( line ) );
        }
        reader.close();
        log.info( "Row names size:" + imageSeriesNames.size() );
        log.info( "Gene names size:" + getAllGenes().size() );
        log.info( "Region names size:" + regionNames.size() );
        log.info( "lines (minus header):" + linecount );
    }

    public Set<String> getAllGenes() {
        return getAllGenes( imageSeriesNames );
    }

    public static Set<String> getAllGenes( Collection<String> rows ) {
        Set<String> result = new HashSet<String>();
        for ( String rowID : rows ) {
            StringTokenizer tokes = new StringTokenizer( rowID, "[" );
            // ugly, uses string formatting
            result.add( tokes.nextToken() );
        }
        return result;
    }

    protected abstract String getRowName( String[] line );

}
