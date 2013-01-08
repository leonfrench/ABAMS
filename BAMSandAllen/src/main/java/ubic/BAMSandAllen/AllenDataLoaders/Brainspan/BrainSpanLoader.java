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
package ubic.BAMSandAllen.AllenDataLoaders.Brainspan;

import java.io.FileReader;
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.Util;
import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import au.com.bytecode.opencsv.CSVReader;

@Deprecated
// replaced by BrainSpanGeneric.java
public class BrainSpanLoader {
    private static Log log = LogFactory.getLog( BrainSpanLoader.class.getName() );

    // values from exon row meta data.csv
    final int gene_symbolPOSITION = 3;
    final int startPOSITION = 5;
    final int endPOSITION = 6;

    // values from col meta data
    final int donor_namePOSITION = 2;
    final int agePOSITION = 3;
    final int genderPOSITION = 4;
    final int structure_idPOSITION = 7;

    String baseFolder = "/home/leon/Downloads/BrainSpan/exons/";
    String prefix;
    private DoubleMatrix<String[], String[]> data;

    public BrainSpanLoader() throws Exception {
        data = new DenseDoubleMatrix<String[], String[]>( 228214, 579 );
        data.setColumnNames( getLines( baseFolder + "exons_columns_metadata.csv" ) );
        data.setRowNames( getLines( baseFolder + "exons_rows_metadata.csv" ) );
        String dataFile = baseFolder + "exons_matrix.csv";
        CSVReader reader = new CSVReader( new FileReader( dataFile ) );
        String line[];
        int row = 0;
        while ( ( line = reader.readNext() ) != null ) {
            // log.info( line.length );
            if ( row % 10000 == 0 ) {
                log.info( row );
                // break;
            }
            // shift over by one because the first column is the Ensembl ID (probably should check that against the row
            // data
            for ( int i = 1; i < line.length; i++ ) {
                data.set( row, i - 1, Double.parseDouble( line[i] ) );
            }
            row++;
        }
    }

    public List<String[]> getLines( String filename ) throws Exception {
        CSVReader reader = new CSVReader( new FileReader( filename ) );
        List<String[]> result = ( List<String[]> ) reader.readAll();
        // eat the header line
        result.remove( 0 );
        reader.close();
        return result;
    }

    public DoubleMatrix<String, String> getGeneData( String gene ) {
        // make a new matrix?
        DoubleMatrix<String, String> result;

        List<String[]> rowsWithGene = new LinkedList<String[]>();
        for ( String[] row : data.getRowNames() ) {
            if ( row[gene_symbolPOSITION].equals( gene ) ) {
                rowsWithGene.add( row );
            }
        }
        log.info( "Rows/exons for gene:" + rowsWithGene.size() );
        result = new DenseDoubleMatrix<String, String>( rowsWithGene.size(), data.columns() );

        for ( String[] row : rowsWithGene ) {
            String rowName = "exon." + row[startPOSITION] + ".to." + row[endPOSITION] + "." + row[gene_symbolPOSITION];
            result.addRowName( rowName );
            for ( int col = 0; col < data.columns(); col++ ) {
                result.set( result.getRowIndexByName( rowName ), col, data.get( data.getRowIndexByName( row ), col ) );
            }
        }

        // colapse column data
        for ( String[] col : data.getColNames() ) {
            String colName = col[agePOSITION] + "." + col[genderPOSITION] + "." + col[structure_idPOSITION] + "."
                    + col[donor_namePOSITION];
            result.addColumnName( colName );
        }
        return result;
    }

    /**
     * @param args
     */
    public static void main( String[] args ) throws Exception {
        // TODO Auto-generated method stub
        BrainSpanLoader loader = new BrainSpanLoader();
        String gene = "SLC20A2";
        DoubleMatrix<String, String> geneMatrix = loader.getGeneData( gene );
        Util.writeRTable( "/home/leon/Downloads/BrainSpan/" + gene + ".txt", geneMatrix );
    }

}
